"""
Reader for GeneHackman LD packs (format version 1).

An LD pack holds precomputed LD (Pearson r) for one or more ancestries, without any
individual-level data.  Layout (see reference_data/pan_ukb_ld/README.md):

    <pack>/MANIFEST.json
    <pack>/<POP>/info.tsv                     KEY<TAB>VALUE (format_version, scale, radius_bp, n_samples, ...)
    <pack>/<POP>/variants_chr<C>.tsv.gz       IDX, CHR, BP37, BP38, EA, OA, EAF, SNP37, SNP38, SEG_LEN, LIFT_STATUS, ...
    <pack>/<POP>/ld_chr<C>.bin                concatenated xz-compressed chunks of int8 values
    <pack>/<POP>/ld_chr<C>.chunks.tsv         CHUNK, ROW_START, ROW_END, OFFSET, NBYTES, NVALUES

Rows are variants in GRCh37 order.  Row i stores r(i, i+1 .. i+SEG_LEN_i) as int8 (r * scale),
i.e. only the upper triangle within the stored radius; the diagonal is implicitly 1 and pairs
beyond the radius are not stored (returned as 0).  r is oriented to the canonical effect allele
(EA, the alphabetically smaller allele, exactly as GeneHackman's standardise_alleles), so it can
be used directly with standardised GWAS z-scores.  SNP37 / SNP38 are GeneHackman SNP keys
(CHR:BP_EA_OA) in each build.

The R twin of this module is R/ld_reference.R; both are tested against the same fixture.
"""

import lzma
import os

import numpy as np
import pandas as pd

FORMAT_NAME = "genehackman-ldpack"
FORMAT_VERSION = 1
DEFAULT_SCALE = 127
VARIANT_COLUMNS = ["IDX", "CHR", "BP37", "BP38", "EA", "OA", "EAF", "SNP37", "SNP38",
                   "SEG_LEN", "LIFT_STATUS", "REF37", "ALT37", "RSID_UKB"]
CHUNK_COLUMNS = ["CHUNK", "ROW_START", "ROW_END", "OFFSET", "NBYTES", "NVALUES"]

PAN_UKB_ANCESTRIES = ["EUR", "AFR", "AMR", "CSA", "EAS", "MID"]
ANCESTRY_ALIASES = {"SAS": "CSA"}
BUILD_KEY_COLUMNS = {"GRCh37": "SNP37", "GRCh38": "SNP38"}
BUILD_BP_COLUMNS = {"GRCh37": "BP37", "GRCh38": "BP38"}


def resolve_ancestry(ancestry):
    """Map a GeneHackman ancestry code onto the pack's population directory name."""
    code = str(ancestry).strip().upper()
    return ANCESTRY_ALIASES.get(code, code)


def key_column(build):
    if build not in BUILD_KEY_COLUMNS:
        raise ValueError(f"Unsupported build for LD pack lookups: {build}. Use GRCh37 or GRCh38.")
    return BUILD_KEY_COLUMNS[build]


def read_info(pack_dir, ancestry):
    """Read <pack>/<POP>/info.tsv into a dict (values left as strings)."""
    pop = resolve_ancestry(ancestry)
    path = os.path.join(pack_dir, pop, "info.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"LD pack info not found for ancestry {ancestry} ({pop}): {path}. "
            "See reference_data/pan_ukb_ld/README.md for how to build or download the pack."
        )
    info = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition("\t")
            info[key] = value
    if int(info.get("format_version", -1)) != FORMAT_VERSION:
        raise ValueError(f"Unsupported LD pack format_version in {path}: {info.get('format_version')}")
    return info


def compute_seg_len(bp37, radius_bp):
    """Number of later variants within radius_bp (inclusive) of each variant; bp37 must be sorted."""
    bp37 = np.asarray(bp37, dtype=np.int64)
    if bp37.size and np.any(np.diff(bp37) < 0):
        raise ValueError("BP37 must be sorted to compute ragged band lengths")
    stop = np.searchsorted(bp37, bp37 + int(radius_bp), side="right") - 1
    return (stop - np.arange(bp37.size)).astype(np.int64)


class LDPack:
    """Random access to one LD pack directory, with per-chromosome caching."""

    def __init__(self, pack_dir):
        self.pack_dir = pack_dir
        self._variants = {}
        self._chunks = {}
        self._info = {}

    def info(self, ancestry):
        pop = resolve_ancestry(ancestry)
        if pop not in self._info:
            self._info[pop] = read_info(self.pack_dir, pop)
        return self._info[pop]

    def scale(self, ancestry):
        return float(self.info(ancestry).get("scale", DEFAULT_SCALE))

    def variants(self, chrom, ancestry):
        pop = resolve_ancestry(ancestry)
        key = (pop, str(chrom))
        if key not in self._variants:
            self.info(pop)
            path = os.path.join(self.pack_dir, pop, f"variants_chr{chrom}.tsv.gz")
            if not os.path.exists(path):
                raise FileNotFoundError(f"LD pack variants not found: {path}")
            text_columns = {"CHR": str, "SNP37": str, "SNP38": str, "EA": str, "OA": str, "LIFT_STATUS": str}
            df = pd.read_csv(path, sep="\t", dtype=text_columns)
            df["IDX"] = df["IDX"].astype(np.int64)
            df["SEG_LEN"] = df["SEG_LEN"].astype(np.int64)
            row_offset = np.zeros(len(df), dtype=np.int64)
            chunks = self.chunks(chrom, pop)
            for _, chunk in chunks.iterrows():
                start, end = int(chunk.ROW_START), int(chunk.ROW_END)
                seg = df["SEG_LEN"].values[start:end]
                row_offset[start:end] = np.concatenate(([0], np.cumsum(seg)[:-1]))
            df["ROW_OFFSET"] = row_offset
            self._variants[key] = df
        return self._variants[key]

    def chunks(self, chrom, ancestry):
        pop = resolve_ancestry(ancestry)
        key = (pop, str(chrom))
        if key not in self._chunks:
            path = os.path.join(self.pack_dir, pop, f"ld_chr{chrom}.chunks.tsv")
            if not os.path.exists(path):
                raise FileNotFoundError(f"LD pack chunk index not found: {path}")
            self._chunks[key] = pd.read_csv(path, sep="\t")
        return self._chunks[key]

    def read_chunk_values(self, chrom, ancestry, chunk_row):
        pop = resolve_ancestry(ancestry)
        path = os.path.join(self.pack_dir, pop, f"ld_chr{chrom}.bin")
        with open(path, "rb") as handle:
            handle.seek(int(chunk_row.OFFSET))
            payload = handle.read(int(chunk_row.NBYTES))
        raw = lzma.decompress(payload)
        values = np.frombuffer(raw, dtype=np.int8)
        if values.size != int(chunk_row.NVALUES):
            raise ValueError(f"Corrupt LD pack chunk {int(chunk_row.CHUNK)} on chr{chrom} ({pop}): "
                             f"expected {int(chunk_row.NVALUES)} values, found {values.size}")
        return values

    def matrix_for_rows(self, chrom, ancestry, rows):
        """Symmetric LD matrix (float64) for pack row indices `rows` (any order)."""
        rows = np.asarray(rows, dtype=np.int64)
        variants = self.variants(chrom, ancestry)
        chunks = self.chunks(chrom, ancestry)
        scale = self.scale(ancestry)
        n = rows.size
        order = np.argsort(rows, kind="stable")
        sorted_rows = rows[order]
        out_sorted = np.eye(n, dtype=np.float64)

        chunk_of_row = np.searchsorted(chunks["ROW_END"].values, sorted_rows, side="right")
        seg_len = variants["SEG_LEN"].values
        row_offset = variants["ROW_OFFSET"].values
        for chunk_index in np.unique(chunk_of_row):
            chunk_row = chunks.iloc[int(chunk_index)]
            values = self.read_chunk_values(chrom, ancestry, chunk_row)
            for pos in np.where(chunk_of_row == chunk_index)[0]:
                i = sorted_rows[pos]
                seg = int(seg_len[i])
                if seg == 0:
                    continue
                stop = np.searchsorted(sorted_rows, i + seg, side="right")
                later = np.arange(pos + 1, stop)
                if later.size == 0:
                    continue
                start = int(row_offset[i])
                r = values[start + (sorted_rows[later] - i - 1)].astype(np.float64) / scale
                out_sorted[pos, later] = r
                out_sorted[later, pos] = r

        inverse = np.empty(n, dtype=np.int64)
        inverse[order] = np.arange(n)
        return out_sorted[np.ix_(inverse, inverse)]

    def get_ld_matrix(self, snps, chrom, ancestry, build):
        """LD matrix for GeneHackman SNP keys in `build`.

        Returns (R, found_snps): R is ordered as found_snps, which keeps the input order of
        `snps` and drops keys that are not in the pack (or failed liftover for GRCh38).
        """
        variants = self.variants(chrom, ancestry)
        column = key_column(build)
        lookup = variants.loc[variants[column].notna(), [column, "IDX"]].drop_duplicates(column)
        index = pd.Series(lookup["IDX"].values, index=lookup[column].values)
        snps = [str(s).upper() for s in snps]
        seen = set()
        found, rows = [], []
        for snp in snps:
            if snp in index.index and snp not in seen:
                seen.add(snp)
                found.append(snp)
                rows.append(int(index[snp]))
        if len(found) == 0:
            return np.zeros((0, 0)), []
        return self.matrix_for_rows(chrom, ancestry, rows), found


def beyond_band_fraction(pack, snps, chrom, ancestry, build):
    """Fraction of off-diagonal pairs among `snps` that lie outside the stored radius."""
    variants = pack.variants(chrom, ancestry)
    column = key_column(build)
    sub = variants[variants[column].isin([str(s).upper() for s in snps])]
    if len(sub) < 2:
        return 0.0
    bp = np.sort(sub["BP37"].values.astype(np.int64))
    radius = int(pack.info(ancestry)["radius_bp"])
    n = bp.size
    within = np.searchsorted(bp, bp + radius, side="right") - np.arange(n) - 1
    total_pairs = n * (n - 1) / 2
    return float(1 - within.sum() / total_pairs)


def make_positive_semidefinite(r_matrix, tolerance=-1e-6):
    """Add (|lambda_min| + 1e-6) to the diagonal when R has a negative eigenvalue below tolerance.

    Returns (R, ridge) where ridge is 0 when no change was needed.
    """
    if r_matrix.shape[0] < 2:
        return r_matrix, 0.0
    min_eig = float(np.linalg.eigvalsh(r_matrix).min())
    if min_eig >= tolerance:
        return r_matrix, 0.0
    ridge = abs(min_eig) + 1e-6
    return r_matrix + ridge * np.eye(r_matrix.shape[0]), ridge
