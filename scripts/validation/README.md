# Harmonisation concordance with MungeSumstats

For the revision 1 response (Reviewer 1 Major 3, Reviewer 2 comment 1): run GeneHackman and MungeSumstats on the
same GWAS and compare rsIDs, effect directions and effect sizes. The comparison logic is `compare_harmonisation()`
in `R/validation.R` (unit-tested in `tests/testthat/test_validation.R`).

This needs the GeneHackman reference data (dbSNP of the output build, liftover chains) and MungeSumstats with its
reference packages (about 6 GB), so it is run on ieup1, not in CI.

## 1. Install MungeSumstats into a user library for the container's R

```bash
SIF=/local-scratch/data/pan_ukb_ld/genehackman_develop.sif   # or the 1.2.0 image
export R_LIBS_USER=/local-scratch/data/genehackman_validation/Rlib
mkdir -p "$R_LIBS_USER"
apptainer exec --env R_LIBS_USER="$R_LIBS_USER" "$SIF" Rscript -e '
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP155.GRCh37", "SNPlocs.Hsapiens.dbSNP155.GRCh38",
                         "BSgenome.Hsapiens.1000genomes.hs37d5", "BSgenome.Hsapiens.NCBI.GRCh38"),
                       update = FALSE, ask = FALSE)'
```

## 2. Run (from `scripts/`, with `.env` pointing at the 1.2.0 reference data)

Three runs, each a few hours and up to ~100 GB RAM for a 10-million-variant file:

```bash
cd ~/GeneHackman/scripts
run() { apptainer exec --env R_LIBS_USER="$R_LIBS_USER" -B /local-scratch "$SIF" \
          Rscript validation/compare_with_mungesumstats.R --threads 8 "$@"; }
OUT=/local-scratch/data/genehackman_validation

# HERMES HFall (GRCh37 input, METAL columns)
run --gwas <HFall file> --columns metal --N <N> --input_build GRCh37 --output_build GRCh38 --output_prefix $OUT/hfall

# Biobank Japan HF (as reformatted for the paper: CHR BP EA OA BETA SE EAF P N; check its build first)
run --gwas hf_japan_GCST90668009.tsv.gz --input_build <GRCh37|GRCh38> --output_build GRCh38 --output_prefix $OUT/bbj

# Reviewer 1's case: the same file identified only by chr:pos:a1:a2, harmonised to GRCh38 with rsIDs
run --gwas hf_japan_GCST90668009.tsv.gz --chrpos_only --input_build <GRCh37|GRCh38> --output_build GRCh38 \
    --output_prefix $OUT/bbj_chrpos
```

## 3. What to report

From each `<prefix>_summary.tsv`: `variants_matched`, `rsid_identical_pct`, `effect_direction_concordant_pct` and
`aligned_beta_correlation`, plus `only_genehackman` / `only_other` (variants dropped by one tool). Look at
`<prefix>_discordant.tsv.gz` before reporting: expected differences are multi-allelic sites (rsID choice),
palindromic variants and variants either tool removes in liftover. The `bbj_chrpos` run should give the same
numbers as `bbj`; if GeneHackman assigns no rsIDs there, the reviewer's bug is not fixed.
