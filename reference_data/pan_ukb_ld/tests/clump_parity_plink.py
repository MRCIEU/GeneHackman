"""
Parity test: GeneHackman's LD-pack clumper against PLINK 1.9 --clump on the same genotypes.

    clump_parity_plink.py --plink plink1.9 --rscript Rscript [--bfile <prefix> --chrom 22] [--out-dir d]

Without --bfile it simulates genotypes (500 people, ~1,500 variants on chr22 with block-like LD);
with --bfile it uses a real panel (e.g. 1000G EUR) restricted to --chrom and --window.  It builds an
LD pack from PLINK's own r (--r, ±RADIUS, counted on A1), simulates a GWAS with several causal
variants, then clumps with both tools under several settings.  Clumps must be identical except where
a pair's r^2 is within int8 quantisation of the threshold.  Two exact checks: an independent
reimplementation of PLINK's rules on PLINK's own r gives PLINK's clumps (rules right), and the R clumper
gives exactly what those rules give on the decoded pack values (implementation right).
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from ldpack_writer import ChromosomeWriter, base_info, write_info  # noqa: E402

RADIUS_KB = 500
SETTINGS = [(0.1, 250), (0.01, 500), (0.5, 500), (0.2, 100)]


def simulate_bfile(prefix, rng, n=500, m=1500, chrom="22"):
    bp = 16_050_000 + np.cumsum(rng.integers(200, 4000, size=m))
    block = np.cumsum(rng.random(m) < 0.04)  # LD blocks of ~25 variants
    latent = np.zeros((2 * n, m))
    for b in np.unique(block):
        cols = np.where(block == b)[0]
        shared = rng.standard_normal((2 * n, 1))
        rho = rng.uniform(0.3, 0.97)
        latent[:, cols] = rho * shared + np.sqrt(1 - rho ** 2) * rng.standard_normal((2 * n, cols.size))
    freq = rng.uniform(0.03, 0.5, size=m)
    haps = latent > np.array([np.quantile(latent[:, j], 1 - freq[j]) for j in range(m)])
    geno = haps[:n].astype(int) + haps[n:].astype(int)  # copies of the "A1" allele
    bases = np.array(list("ACGT"))
    a1 = rng.choice(bases, size=m)
    a2 = np.array([rng.choice(bases[bases != x]) for x in a1])
    with open(prefix + ".map", "w") as handle:
        for j in range(m):
            handle.write(f"{chrom}\trs{j + 1}\t0\t{bp[j]}\n")
    with open(prefix + ".ped", "w") as handle:
        for i in range(n):
            alleles = []
            for j in range(m):
                g = geno[i, j]
                alleles += [a1[j], a1[j]] if g == 2 else ([a1[j], a2[j]] if g == 1 else [a2[j], a2[j]])
            handle.write(f"F{i} I{i} 0 0 0 -9 " + " ".join(alleles) + "\n")
    return a1


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plink", default="plink1.9")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--bfile", default=None)
    parser.add_argument("--chrom", default="22")
    parser.add_argument("--window", default=None, help="START-END to restrict a real panel")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--ld", choices=["haplotype", "genotype"], default="haplotype",
                        help="r used to build the pack (haplotype = what PLINK --clump uses)")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    work = args.out_dir or tempfile.mkdtemp(prefix="clump_parity_")
    os.makedirs(work, exist_ok=True)
    bfile = os.path.join(work, "panel")

    if args.bfile:
        extract = ["--bfile", args.bfile, "--chr", args.chrom]
        if args.window:
            start, end = args.window.split("-")
            extract += ["--from-bp", start, "--to-bp", end]
        run([args.plink] + extract + ["--snps-only", "just-acgt", "--biallelic-only", "strict", "--maf", "0.01",
                                      "--keep-allele-order", "--make-bed", "--out", bfile])
    else:
        simulate_bfile(bfile, rng, chrom=args.chrom)
        run([args.plink, "--file", bfile, "--keep-allele-order", "--make-bed", "--out", bfile])

    bim = pd.read_csv(bfile + ".bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                      dtype={"CHR": str, "A1": str, "A2": str})
    bim = bim.drop_duplicates("BP", keep=False).reset_index(drop=True)  # one variant per position, for simple keys
    bim["SNP"].to_csv(os.path.join(work, "keep.txt"), index=False, header=False)
    run([args.plink, "--bfile", bfile, "--extract", os.path.join(work, "keep.txt"), "--keep-allele-order",
         "--make-bed", "--out", bfile + "_u"])
    bfile = bfile + "_u"
    bim = pd.read_csv(bfile + ".bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                      dtype={"CHR": str, "A1": str, "A2": str})
    m = len(bim)

    # r from PLINK (counted on A1) within the radius, and an LD pack built from it.  PLINK's --clump uses
    # haplotype-based r^2 (EM, as reported by --r dprime); "--ld genotype" uses genotype correlations
    # (PLINK --r, which is what Pan-UKB's matrices are) to show how much that choice alone changes clumps.
    ld_flags = ["--r", "dprime"] if args.ld == "haplotype" else ["--r"]
    run([args.plink, "--bfile", bfile] + ld_flags
        + ["--ld-window-kb", str(RADIUS_KB), "--ld-window", "999999", "--ld-window-r2", "0", "--keep-allele-order",
           "--out", os.path.join(work, "r")])
    pairs = pd.read_csv(os.path.join(work, "r.ld"), sep=r"\s+")
    position = {s: i for i, s in enumerate(bim["SNP"])}
    r_a1 = np.eye(m)
    ia, ib = pairs["SNP_A"].map(position).values, pairs["SNP_B"].map(position).values
    r_a1[ia, ib] = np.nan_to_num(pairs["R"].values)
    r_a1[ib, ia] = np.nan_to_num(pairs["R"].values)
    ea = np.where(bim["A1"] < bim["A2"], bim["A1"], bim["A2"])
    oa = np.where(bim["A1"] < bim["A2"], bim["A2"], bim["A1"])
    variants = pd.DataFrame({
        "CHR": args.chrom, "BP37": bim["BP"].values, "EA": ea, "OA": oa, "EAF": 0.5,
        "SNP37": [f"{args.chrom}:{b}_{e}_{o}" for b, e, o in zip(bim["BP"], ea, oa)],
        "SNP38": None, "BP38": None, "LIFT_STATUS": "unmapped", "REF37": bim["A2"].values, "ALT37": bim["A1"].values,
        "SIGN": np.where(ea == bim["A1"].values, 1, -1), "IDX_GLOBAL": np.arange(m)})
    pack = os.path.join(work, "pack")
    writer = ChromosomeWriter(os.path.join(pack, "EUR"), os.path.join(work, "chunks"), args.chrom, variants,
                              RADIUS_KB * 1000, chunk_rows=256, deadzone=0.0)
    for k in writer.pending_chunks():
        r0, r1 = writer.bounds[k]
        writer.write_chunk(k, r_a1[r0:r1, r0:writer.column_end(k)])
    writer.finalise()
    write_info(os.path.join(pack, "EUR"), base_info(500, RADIUS_KB * 1000, 0.01, 0.0, "PLINK parity test", 256))

    # GWAS: z = R b + noise, with a few causal variants
    sign = variants["SIGN"].values
    r_ea = r_a1 * np.outer(sign, sign)
    beta = np.zeros(m)
    beta[rng.choice(m, size=6, replace=False)] = rng.choice([-1, 1], size=6) * rng.uniform(8, 14, size=6)
    z = np.clip(r_ea, -1, 1) @ beta + rng.standard_normal(m)
    from scipy.stats import norm
    p = 2 * norm.sf(np.abs(z))
    gwas = pd.DataFrame({"SNP": variants["SNP37"], "CHR": args.chrom, "BP": bim["BP"], "EA": ea, "OA": oa,
                         "BETA": z / 100, "SE": 0.01, "P": p, "RSID": bim["SNP"]})
    gwas_file = os.path.join(work, "gwas.tsv")
    gwas.to_csv(gwas_file, sep="\t", index=False)
    gwas[["RSID", "P"]].rename(columns={"RSID": "SNP"}).to_csv(os.path.join(work, "assoc.txt"), sep="\t", index=False)
    print(f"{m} variants, {int((p < 5e-8).sum())} with P < 5e-8, {int((p < 0.01).sum())} with P < 0.01")

    # Decoded pack values: what the R clumper sees after int8 quantisation
    sys.path.insert(0, os.path.join(REPO, "scripts"))
    from ld_reference import LDPack
    r_pack = LDPack(pack).matrix_for_rows(args.chrom, "EUR", np.arange(m)) * np.outer(sign, sign)  # back to A1

    def members(frame):
        out_map = {}
        for row in frame.itertuples():
            out_map[row.SNP] = set() if row.SP2 == "NONE" else {s_.split("(")[0] for s_ in row.SP2.split(",")}
        return out_map

    def greedy(r_matrix, r2, kb, p1=5e-8, p2=0.01):
        """Independent reimplementation of PLINK 1.9 --clump, on a given r matrix."""
        pool = np.where(p <= p2)[0]
        order = pool[np.lexsort((bim["BP"].values[pool], p[pool]))]
        clumped = np.zeros(m, dtype=bool)
        clumps = {}
        for i in order:
            if p[i] > p1 or clumped[i]:
                continue
            clumped[i] = True
            near = pool[(~clumped[pool]) & (np.abs(bim["BP"].values[pool] - bim["BP"].values[i]) <= kb * 1000)]
            joined = near[r_matrix[i, near] ** 2 >= r2]
            clumped[joined] = True
            clumps[bim["SNP"][i]] = set(bim["SNP"][joined])
        return clumps

    failures = 0
    for r2, kb in SETTINGS:
        clump_args = f"--clump-p1 5e-8 --clump-p2 0.01 --clump-r2 {r2} --clump-kb {kb}"
        out = os.path.join(work, f"plink_{r2}_{kb}")
        run([args.plink, "--bfile", bfile, "--clump", os.path.join(work, "assoc.txt")] + clump_args.split()
            + ["--out", out])
        plink = members(pd.read_csv(out + ".clumped", sep=r"\s+")) if os.path.exists(out + ".clumped") else {}
        ours_file = os.path.join(work, f"ours_{r2}_{kb}.clumped")
        code = ("suppressMessages(pkgload::load_all('{repo}', quiet = TRUE));"
                "invisible(clump_gwas_with_ld_pack('{gwas}', 'EUR', 'GRCh37', ' {args}', '{out}', '{pack}'))")
        r_code = code.format(repo=REPO, gwas=gwas_file, args=clump_args, out=ours_file, pack=pack)
        subprocess.run([args.rscript, "-e", r_code], check=True, stderr=subprocess.DEVNULL)
        ours = members(pd.read_csv(ours_file, sep=r"\s+"))

        rules_match_plink = greedy(r_a1, r2, kb) == plink         # same rules, same (float) r
        implementation_ok = greedy(r_pack, r2, kb) == ours        # R clumper = rules on decoded pack r
        same = sum(plink.get(k) == ours.get(k) for k in plink)
        if args.ld == "haplotype":
            ok = rules_match_plink and implementation_ok
        else:
            ok = implementation_ok
        print(f"{'OK  ' if ok else 'FAIL'} r2 {r2}, kb {kb}: rules reproduce PLINK on its own r: {rules_match_plink}; "
              f"R clumper reproduces the rules on the pack: {implementation_ok}; "
              f"clumps identical to PLINK: {same} of {len(plink)} (ours {len(ours)})")
        failures += not ok
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
