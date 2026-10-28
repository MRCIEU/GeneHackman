"""
Phase 0 check: compare a pack chromosome with the 1000 Genomes panel GeneHackman already uses.

    61_compare_1000g.py --pack <pack> --pop EUR --chrom 22 \
        --bfile <genomic_data>/1000genomes/b37_dbsnp156/EUR [--frq <...>/EUR.frq] \
        --window 20000000-25000000 [--plink plink1.9] [--out compare_EUR_chr22.json]

Checks that (1) the pack's EAF is the frequency of the same allele as in 1000G (no allele swaps:
the pack AF orientation is undocumented upstream), and (2) r has the same sign as PLINK's r on the
1000G genotypes for clearly correlated pairs (|r| > 0.3), i.e. the pack's orientation to the
canonical effect allele is right.  Variants are matched on the GeneHackman GRCh37 key, built from
the bim alleles (EA = alphabetically smaller allele).  Works anywhere PLINK 1.9 and the 1000G files
are available (it does not need Hail or the source matrices).
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from ld_reference import LDPack  # noqa: E402


def canonical_key(chrom, bp, a1, a2):
    a1, a2 = a1.upper(), a2.upper()
    ea, oa = (a1, a2) if a1 < a2 else (a2, a1)
    return f"{chrom}:{bp}_{ea}_{oa}", ea


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pack", required=True)
    parser.add_argument("--pop", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--frq", default=None)
    parser.add_argument("--window", default=None, help="START-END (GRCh37) for the r comparison")
    parser.add_argument("--plink", default="plink1.9")
    parser.add_argument("--min-abs-r", type=float, default=0.3)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    pack = LDPack(args.pack)
    variants = pack.variants(args.chrom, args.pop)
    radius = int(pack.info(args.pop)["radius_bp"])
    bim = pd.read_csv(args.bfile + ".bim", sep=r"\s+", header=None, names=["CHR", "RSID", "CM", "BP", "A1", "A2"],
                      dtype={"CHR": str, "A1": str, "A2": str})
    bim = bim[bim["CHR"] == str(args.chrom)].copy()
    keys = [canonical_key(c, b, a1, a2) for c, b, a1, a2 in zip(bim["CHR"], bim["BP"], bim["A1"], bim["A2"])]
    bim["KEY"] = [k for k, _ in keys]
    bim["EA"] = [e for _, e in keys]
    bim["A1_IS_EA"] = bim["A1"].str.upper() == bim["EA"]
    merged = bim.merge(variants[["SNP37", "IDX", "BP37", "EAF"]], left_on="KEY", right_on="SNP37")
    result = {"pop": args.pop, "chrom": args.chrom, "n_pack": int(len(variants)), "n_1000g": int(len(bim)),
              "n_overlap": int(len(merged))}

    if args.frq and os.path.exists(args.frq):
        frq = pd.read_csv(args.frq, sep=r"\s+", dtype={"A1": str, "A2": str})
        frq = frq[frq["SNP"].isin(merged["RSID"])].drop_duplicates("SNP").set_index("SNP")
        m = merged[merged["RSID"].isin(frq.index)].copy()
        a1 = frq.loc[m["RSID"], "A1"].str.upper().values
        maf = frq.loc[m["RSID"], "MAF"].values.astype(float)
        m["EAF_1000G"] = np.where(a1 == m["EA"].values, maf, 1 - maf)
        diff = m["EAF"] - m["EAF_1000G"]
        swapped = np.abs(m["EAF"] - (1 - m["EAF_1000G"])) < np.abs(diff)
        result["eaf"] = {"n": int(len(m)), "pearson": float(np.corrcoef(m["EAF"], m["EAF_1000G"])[0, 1]),
                         "mean_abs_diff": float(np.abs(diff).mean()),
                         "n_abs_diff_gt_0.2": int((np.abs(diff) > 0.2).sum()),
                         "share_closer_if_swapped": float(swapped.mean())}

    window = merged
    if args.window:
        start, end = (int(x) for x in args.window.split("-"))
        window = merged[(merged["BP"] >= start) & (merged["BP"] <= end)]
    window = window.sort_values("BP").drop_duplicates("RSID").drop_duplicates("KEY")
    with tempfile.TemporaryDirectory() as tmp:
        snp_file = os.path.join(tmp, "snps.txt")
        window["RSID"].to_csv(snp_file, index=False, header=False)
        out = os.path.join(tmp, "ld")
        subprocess.run([args.plink, "--bfile", args.bfile, "--chr", str(args.chrom), "--extract", snp_file,
                        "--r", "--ld-window-kb", str(radius // 1000), "--ld-window", "999999",
                        "--ld-window-r2", "0", "--keep-allele-order", "--out", out],
                       check=True, stdout=subprocess.DEVNULL)
        pairs = pd.read_csv(out + ".ld", sep=r"\s+")
    lookup = window.set_index("RSID")
    pairs = pairs[pairs["SNP_A"].isin(lookup.index) & pairs["SNP_B"].isin(lookup.index)]
    sign_a = np.where(lookup.loc[pairs["SNP_A"], "A1_IS_EA"].values, 1, -1)
    sign_b = np.where(lookup.loc[pairs["SNP_B"], "A1_IS_EA"].values, 1, -1)
    r_1000g = pairs["R"].values * sign_a * sign_b
    keys = list(dict.fromkeys(lookup["KEY"].values))
    r_pack_matrix, found = pack.get_ld_matrix(keys, args.chrom, args.pop, "GRCh37")
    position = {k: i for i, k in enumerate(found)}
    key_a = lookup.loc[pairs["SNP_A"], "KEY"].values
    key_b = lookup.loc[pairs["SNP_B"], "KEY"].values
    r_pack = np.array([r_pack_matrix[position[a], position[b]] for a, b in zip(key_a, key_b)])
    strong = np.abs(r_1000g) > args.min_abs_r
    result["r"] = {"n_pairs": int(len(r_pack)), "pearson": float(np.corrcoef(r_pack, r_1000g)[0, 1]),
                   "n_strong_pairs": int(strong.sum()),
                   "sign_concordance_strong": float((np.sign(r_pack[strong]) == np.sign(r_1000g[strong])).mean())
                   if strong.any() else None,
                   "mean_abs_diff": float(np.abs(r_pack - r_1000g).mean())}
    text = json.dumps(result, indent=1)
    print(text)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
