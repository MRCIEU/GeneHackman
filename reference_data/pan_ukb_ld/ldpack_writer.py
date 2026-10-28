"""
Writer for GeneHackman LD packs (format version 1); the reader is scripts/ld_reference.py.

The writer is source-agnostic: callers hand it, chunk by chunk, a dense float rectangle of r
(rows r0..r1-1, columns r0..c_end-1, oriented to the ALT allele as stored by Pan-UKB/gnomAD/PLINK
counted allele) and it applies the per-variant SIGN (EA orientation), the dead-zone and int8
quantisation, xz-compresses, and checkpoints each chunk to disk so a multi-day build can resume.
"""

import json
import lzma
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from ld_reference import (CHUNK_COLUMNS, DEFAULT_SCALE, FORMAT_NAME, FORMAT_VERSION,  # noqa: E402
                          VARIANT_COLUMNS, compute_seg_len)

XZ_PRESET = 6


def chunk_bounds(n_rows, chunk_rows):
    return [(start, min(start + chunk_rows, n_rows)) for start in range(0, n_rows, chunk_rows)]


def chunk_column_end(r0, r1, seg_len):
    """Exclusive end column needed to cover the ragged band of rows r0..r1-1."""
    rows = np.arange(r0, r1)
    return int(max(r1, (rows + seg_len[r0:r1]).max() + 1))


def quantise(r, scale=DEFAULT_SCALE, deadzone=0.0):
    r = np.nan_to_num(np.asarray(r, dtype=np.float64), nan=0.0)
    q = np.clip(np.rint(r * scale), -scale, scale)
    if deadzone > 0:
        q[np.abs(r) < deadzone] = 0
    return q.astype(np.int8)


def encode_chunk(rect, r0, r1, seg_len, sign, scale=DEFAULT_SCALE, deadzone=0.0):
    """Encode rows r0..r1-1 of the ragged band from `rect` (rows r0..r1-1 x cols r0..)."""
    pieces = []
    n_nan = 0
    for i in range(r0, r1):
        seg = int(seg_len[i])
        if seg == 0:
            continue
        local_row = i - r0
        cols = np.arange(i + 1, i + seg + 1)
        r = np.asarray(rect[local_row, cols - r0], dtype=np.float64)
        n_nan += int(np.isnan(r).sum())
        r = r * sign[i] * sign[cols]
        pieces.append(quantise(r, scale, deadzone))
    values = np.concatenate(pieces) if pieces else np.zeros(0, dtype=np.int8)
    return values, n_nan


class ChromosomeWriter:
    """Writes one chromosome of one population, checkpointing each chunk under work_dir."""

    def __init__(self, pack_pop_dir, work_dir, chrom, variants, radius_bp, chunk_rows=2048,
                 scale=DEFAULT_SCALE, deadzone=0.0):
        required = {"CHR", "BP37", "EA", "OA", "EAF", "SNP37", "SNP38", "BP38", "LIFT_STATUS",
                    "REF37", "ALT37", "SIGN"}
        missing = required - set(variants.columns)
        if missing:
            raise ValueError(f"variants table missing columns: {sorted(missing)}")
        # Row order is GRCh37 position; ties (several alleles at one position) follow the source
        # matrix index when known, so the rows requested from the source stay strictly increasing.
        tie_breaker = "IDX_GLOBAL" if "IDX_GLOBAL" in variants.columns else "SNP37"
        self.variants = variants.sort_values(["BP37", tie_breaker], kind="stable").reset_index(drop=True)
        if "IDX_GLOBAL" in self.variants.columns and np.any(np.diff(self.variants["IDX_GLOBAL"].values) <= 0):
            raise ValueError("Source matrix indices are not increasing with GRCh37 position")
        self.pack_pop_dir = pack_pop_dir
        self.chunk_dir = os.path.join(work_dir, f"chunks_chr{chrom}")
        self.chrom = str(chrom)
        self.radius_bp = int(radius_bp)
        self.chunk_rows = int(chunk_rows)
        self.scale = int(scale)
        self.deadzone = float(deadzone)
        self.seg_len = compute_seg_len(self.variants["BP37"].values, self.radius_bp)
        self.sign = self.variants["SIGN"].values.astype(np.int64)
        if not np.all(np.isin(self.sign, [-1, 1])):
            raise ValueError("SIGN must be +1 or -1 for every variant")
        self.bounds = chunk_bounds(len(self.variants), self.chunk_rows)
        os.makedirs(self.chunk_dir, exist_ok=True)
        os.makedirs(self.pack_pop_dir, exist_ok=True)
        self._write_settings()

    def _settings(self):
        return {"chrom": self.chrom, "n_rows": int(len(self.variants)), "radius_bp": self.radius_bp,
                "chunk_rows": self.chunk_rows, "scale": self.scale, "deadzone": self.deadzone}

    def _write_settings(self):
        path = os.path.join(self.chunk_dir, "settings.json")
        settings = self._settings()
        if os.path.exists(path):
            with open(path, encoding="utf-8") as handle:
                previous = json.load(handle)
            if previous != settings:
                raise ValueError(
                    f"Existing checkpoints in {self.chunk_dir} were written with different settings "
                    f"({previous}) than requested ({settings}). Remove that directory to rebuild.")
        else:
            _atomic_write_text(path, json.dumps(settings, indent=1))

    def chunk_path(self, k):
        return os.path.join(self.chunk_dir, f"chunk_{k:06d}.xz")

    def is_chunk_done(self, k):
        return os.path.exists(self.chunk_path(k) + ".json")

    def pending_chunks(self):
        return [k for k in range(len(self.bounds)) if not self.is_chunk_done(k)]

    def column_end(self, k):
        r0, r1 = self.bounds[k]
        return chunk_column_end(r0, r1, self.seg_len)

    def write_chunk(self, k, rect):
        r0, r1 = self.bounds[k]
        c_end = self.column_end(k)
        if rect.shape[0] != r1 - r0 or rect.shape[1] < c_end - r0:
            raise ValueError(f"chunk {k}: rectangle {rect.shape} does not cover rows {r0}:{r1}, cols {r0}:{c_end}")
        values, n_nan = encode_chunk(rect, r0, r1, self.seg_len, self.sign, self.scale, self.deadzone)
        expected = int(self.seg_len[r0:r1].sum())
        if values.size != expected:
            raise AssertionError(f"chunk {k}: encoded {values.size} values, expected {expected}")
        payload = lzma.compress(values.tobytes(), preset=XZ_PRESET)
        _atomic_write_bytes(self.chunk_path(k), payload)
        meta = {"chunk": k, "row_start": r0, "row_end": r1, "nbytes": len(payload),
                "nvalues": int(values.size), "n_nan": n_nan,
                "n_zero": int((values == 0).sum())}
        _atomic_write_text(self.chunk_path(k) + ".json", json.dumps(meta))
        return meta

    def finalise(self):
        """Concatenate checkpointed chunks into ld_chr<C>.bin / .chunks.tsv and write the variants table."""
        pending = self.pending_chunks()
        if pending:
            raise RuntimeError(f"chr{self.chrom}: {len(pending)} chunks still pending, cannot finalise")
        bin_path = os.path.join(self.pack_pop_dir, f"ld_chr{self.chrom}.bin")
        rows = []
        offset = 0
        tmp_bin = bin_path + ".tmp"
        with open(tmp_bin, "wb") as out:
            for k in range(len(self.bounds)):
                with open(self.chunk_path(k) + ".json", encoding="utf-8") as handle:
                    meta = json.load(handle)
                with open(self.chunk_path(k), "rb") as handle:
                    payload = handle.read()
                if len(payload) != meta["nbytes"]:
                    raise RuntimeError(f"chunk {k} size mismatch on disk")
                out.write(payload)
                rows.append([k, meta["row_start"], meta["row_end"], offset, meta["nbytes"], meta["nvalues"]])
                offset += meta["nbytes"]
        os.replace(tmp_bin, bin_path)
        chunks = pd.DataFrame(rows, columns=CHUNK_COLUMNS)
        _atomic_write_df(chunks, os.path.join(self.pack_pop_dir, f"ld_chr{self.chrom}.chunks.tsv"))

        out = self.variants.copy()
        out["IDX"] = np.arange(len(out))
        out["SEG_LEN"] = self.seg_len
        for column in VARIANT_COLUMNS:
            if column not in out.columns:
                out[column] = pd.NA
        for column in ("IDX", "BP37", "BP38", "SEG_LEN"):
            out[column] = pd.to_numeric(out[column]).round().astype("Int64")
        _atomic_write_df(out[VARIANT_COLUMNS],
                         os.path.join(self.pack_pop_dir, f"variants_chr{self.chrom}.tsv.gz"),
                         compression="gzip")
        totals = {"chrom": self.chrom, "n_variants": int(len(out)), "n_chunks": len(rows),
                  "bytes": int(offset), "n_values": int(self.seg_len.sum())}
        return totals


def write_info(pack_pop_dir, info):
    lines = [f"{key}\t{value}" for key, value in info.items()]
    _atomic_write_text(os.path.join(pack_pop_dir, "info.tsv"), "\n".join(lines) + "\n")


def base_info(n_samples, radius_bp, maf_min, deadzone, source, chunk_rows, scale=DEFAULT_SCALE):
    return {"format": FORMAT_NAME, "format_version": FORMAT_VERSION, "scale": scale,
            "compression": "xz", "radius_bp": int(radius_bp), "maf_min": maf_min,
            "deadzone": deadzone, "n_samples": n_samples, "chunk_rows": chunk_rows,
            "index_build": "GRCh37", "key_builds": "GRCh37,GRCh38", "source": source}


def _atomic_write_bytes(path, payload):
    directory = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=".tmp_")
    with os.fdopen(fd, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, path)


def _atomic_write_text(path, text):
    _atomic_write_bytes(path, text.encode("utf-8"))


def _atomic_write_df(df, path, compression=None):
    directory = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=".tmp_")
    os.close(fd)
    df.to_csv(tmp, sep="\t", index=False, compression=compression)
    os.replace(tmp, path)
