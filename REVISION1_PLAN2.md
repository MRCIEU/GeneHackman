# GeneHackman revision 1, plan 2: run everything on ieup1

**6 November 2026.** Follows `REVISION1_PLAN.md` (reference data, parts A1–B2, all done) and `REVISION1_GAP_ANALYSIS.md` (what each reviewer point needs). This plan covers what still has to run on ieup1:

- the container with the new code;
- the 1.2.0 release;
- the heart-failure (HF) re-analysis that fills the yellow placeholders in the tracked manuscript and response letter;
- the LDSC benchmark;
- the MungeSumstats comparison, including Reviewer 1's `chr:pos:a1:a2` case.

Each part can be started on its own ("start part C4"), in this order:

| Part | What | Needs | Rough time |
|---|---|---|---|
| **C0** | Get the code onto ieup1 | Branch pushed, or a bundle | 10 min |
| **C1** | Container 1.2.0: EnsDb v86, version checks | C0 | 1–2 h (build) |
| **C2** | Unit tests, Python tests and dry runs on ieup1 | C0, C1 | 30 min |
| **C3** | HF data in place; fill in case and control counts | C0 | 30 min plus reading |
| **C4** | HF re-analysis: 14 pipeline runs | C1–C3 | 1–2 days of machine time |
| **C5** | LDSC benchmark (R1-M1) | C4 run 1 (standardised GWAS) | 1 h |
| **C6** | MungeSumstats concordance (R1-M3, R2-1), including Reviewer 1's case | C1, C3 | 1 day of machine time |
| **C7** | Release: upload 1.2.0, end-to-end tests, `testing_complete.txt` | C1, C2 | 1 day |
| **C8** | Send the results back for the documents | C4–C6 | – |

Notes:
- **Container binds.** The pipelines bind the repository's `R/` and `scripts/` into the container (`snakemake/profiles/local/config.yaml`). So the container only needs rebuilding for packages; code changes take effect immediately.
- **Order.** C5 and C6 can run while C4 is running. C7 can wait until after the re-analysis, but the manuscript should cite the 1.2.0 image (C1).

---

## C0. Get the code onto ieup1

All the work is on `feature/1.1.1-reviewer-feedback`. You pushed it up to `f606718`; the commits after that (the LD-mismatch report work, this plan and the HF configs) are on the laptop and not pushed yet.

```bash
# laptop, after reviewing
cd ~/Documents/Projects/GeneHackman && git push origin feature/1.1.1-reviewer-feedback
# ieup1
cd /home/wt23152/GeneHackman
git status                                   # the GRCh38 config fixes from REVISION1_PLAN §A2.1 should already be in 0e4b2e1
git fetch origin && git checkout feature/1.1.1-reviewer-feedback && git pull --ff-only
git log --oneline -1                         # expect the "REVISION1_PLAN2" commit (or later)
```

If ieup1 has uncommitted changes (e.g. `reference_data/grch38/local_config.sh`, which is git-ignored), `git status` shows them. Commit or stash anything else before pulling.

## C1. Container 1.2.0

**Why.** LocusZoom in GRCh38 needs `EnsDb.Hsapiens.v86`, which is already in `docker/requirements.R`. Bumping the version also gives the paper an immutable image tag to cite (audit B32).

1. **Version.** Set `Version: 1.2.0` in `DESCRIPTION` and commit.
2. **Build and push** from a machine with Docker. `--platform linux/amd64` is needed on the Mac.
   ```bash
   docker build --platform linux/amd64 -f docker/Dockerfile -t mrcieu/genehackman:1.2.0 .
   docker push mrcieu/genehackman:1.2.0
   docker inspect --format '{{index .RepoDigests 0}}' mrcieu/genehackman:1.2.0   # the digest to cite in the paper
   ```
3. **SIF on ieup1.** `run_pipeline.sh` builds the SIF itself if it's missing and the directory is writable. To build it by hand:
   ```bash
   apptainer build /home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif docker://mrcieu/genehackman:1.2.0
   ```
4. **Check the packages the new code relies on:**
   ```bash
   SIF=/home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif
   apptainer exec $SIF Rscript -e 'for (p in c("susieR","coloc","EnsDb.Hsapiens.v86","EnsDb.Hsapiens.v75","data.table","knitr","prettydoc")) cat(p, as.character(packageVersion(p)), "\n"); stopifnot(exists("kriging_rss", asNamespace("susieR")), exists("estimate_s_rss", asNamespace("susieR")))'
   apptainer exec $SIF python -c "import MultiSuSiE, scipy, numpy, pandas; print('python ok', scipy.__version__)"
   apptainer exec $SIF bash -c 'which plink1.9 liftOver; ls /home/ldsc/ldsc.py'
   ```
   - `kriging_rss` and `estimate_s_rss` need susieR ≥ 0.12; the tests here used 0.14.2.
   - The MultiSuSiE LD-mismatch check needs scipy, which is in `docker/requirements.txt`.
5. **`.env`:** `DOCKER_VERSION=1.2.0`, or unset it now that `DESCRIPTION` says 1.2.0.

**Done when** all three checks print without errors.

## C2. Tests and dry runs on ieup1

```bash
cd /home/wt23152/GeneHackman
SIF=/home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif
apptainer exec -B /local-scratch $SIF Rscript -e 'pkgload::load_all("."); r <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE); df <- as.data.frame(r); cat("tests:", nrow(df), "failed:", sum(df$failed > 0), "errors:", sum(df$error), "\n")'
apptainer exec $SIF python -m unittest discover tests/python
```

**Expected results:**
- **R tests:** 143 test blocks, 0 failed. Some skip when they don't find the Pan-UKB pack at the default path. Point `options(GeneHackman.ld_pack_dir=...)` at `pipeline_data/genomic_data/pan_ukb_ld` to run those too.
- **Python tests:** 20 tests, OK. This includes `TestLdMismatch`, which compares the Python LD-mismatch check with susieR's output.

Then dry-run each HF config (they're written in C3) with `--dry-run`, for example `./run_pipeline.sh snakemake/finemap.smk paper/revision1_hf/finemap.yaml --dry-run`. Any TODO count still at 0 stops the dry run with a message, which is intended.

## C3. HF data and configuration

### C3.1 Input files

The submitted analysis read the files from BluePebble (`/user/work/wt23152/020/`). The configs in `paper/revision1_hf/` use the relative path `hf_data/…`, so the supplementary YAMLs have no local paths (Reviewer 1 m20). Copy the files to ieup1 and link the folder into the repository root (`hf_data` is git-ignored):

```bash
HF_DIR=/local-scratch/data/genehackman_hf            # any folder with ~10 GB free
mkdir -p $HF_DIR && rsync -av bp1:/user/work/wt23152/020/{FORMAT-METAL_Pheno1_EUR.tsv.gz,FORMAT-METAL_Pheno3_EUR.tsv.gz,FORMAT-METAL_Pheno4_EUR.tsv.gz,hf_progression_allcause_death.autosomes.gz} $HF_DIR/
rsync -av <laptop>:~/Documents/Projects/GeneHackman/hf_japan_GCST90668009.tsv.gz $HF_DIR/   # BBJ; also in the laptop repo folder
cd /home/wt23152/GeneHackman && ln -sfn $HF_DIR hf_data
```

### C3.2 Check the headers and builds

```bash
for f in hf_data/*.gz; do echo "== $f"; zcat $f | head -2; done
```

- **HERMES** (`FORMAT-METAL_Pheno*`): the configs map `rsID chr pos_b37 A1 A2 A1_freq A1_beta se pval` (as in the submitted runs; GRCh37).
  - If the header also has per-variant case and control counts, add `N_CASES: <column>` and `N_CONTROLS: <column>` to the `hermes_columns` block in each config. GeneHackman then computes a per-variant `N_EFF`, which is better for a meta-analysis.
  - If it has a per-variant effective N, map it as `N_EFF`. Either takes precedence over the study totals.
- **BBJ** (`hf_japan_GCST90668009.tsv.gz`): columns `CHR BP EA OA BETA SE EAF P N`, but the build isn't recorded.
  - Check a known variant against dbSNP, e.g. `zcat hf_data/hf_japan_GCST90668009.tsv.gz | awk '$1==1 && $2==591978'` and look up its rsID in both builds.
  - Then set `build:` in `finemap_multi_ancestry*.yaml`. GeneHackman's build check (`liftover.check_build`, on by default) also stops a run if the declared build is wrong.
  - The YAML's N (203,234 in the submitted run) differs from the file's N (206,495). The file's N now takes precedence. Reconcile it for the dataset table: Enzan et al. 2025 report 16,251 cases and 197,577 controls.
- **Progression GWAS** (`hf_progression_allcause_death.autosomes.gz`): check the columns (`rsid chr bp ea oa eaf beta se p`) and the build.

### C3.3 Fill in the TODOs

`grep -n TODO paper/revision1_hf/*.yaml` lists them. A count left at 0 stops the pipeline with a message. Prevalences are marked "decide".

| Value | Configs | Source |
|---|---|---|
| HFall cases and controls | all except `disease_progression_*` (which also use HFall) | Henry et al. 2025, Table 1 / Supplementary Table 1 |
| ni-HFrEF cases and controls | `compare_gwases`, `coloc*` | Henry et al. 2025 |
| ni-HFpEF cases (3,841 in the draft letter, check) and controls | `compare_gwases`, `coloc*` | Henry et al. 2025 |
| HF population prevalence (`prevalence`) | HERMES phenotypes | Decide, e.g. 0.02 for HF and 0.01 for each subtype, and state it in the methods. It changes only the liability-scale h². |
| Deaths and survivors | `disease_progression_*` | Sunderland et al. 2025 |
| BBJ and progression `build:` | as above | C3.2 |

The same counts appear in several files. Edit them with one `sed` per value, for example:

```bash
sed -i 's/N_case: 0  # TODO: HFall cases.*/N_case: <n>  # Henry et al. 2025/' paper/revision1_hf/*.yaml
```

Commit the filled-in configs: they become the supplementary configuration files.

## C4. HF re-analysis

### C4.1 One project folder per variant

Outputs go to `PROJECT_DIR/{data,results}`. Runs that differ only in a setting (L = 10, Pan-UKB LD, CWLS) get their own project. Within a project, standardisation and clumping are shared between pipelines.

```bash
cd /home/wt23152/GeneHackman
OUT=/local-scratch/data/genehackman_hf_runs; mkdir -p $OUT
cp .env .env.backup
run() {   # run <project> <pipeline> <config>: sets PROJECT_DIR, times the run, keeps the log
  local project=$1 smk=$2 yaml=$3; mkdir -p $OUT/$project
  sed -i "s|^PROJECT_DIR=.*|PROJECT_DIR=$OUT/$project|" .env
  /usr/bin/time -v ./run_pipeline.sh snakemake/$smk.smk paper/revision1_hf/$yaml.yaml \
    > $OUT/$project/${smk}_${yaml}.log 2>&1
  echo "$project $smk $yaml exit=$? $(grep -E 'Elapsed|Maximum resident' $OUT/$project/${smk}_${yaml}.log | tr '\n' ' ')" | tee -a $OUT/timings.txt
}
reuse() { mkdir -p $OUT/$2 && cp -al $OUT/$1/data $OUT/$2/; }   # share standardised GWAS and clumps (hard links)
```

- Run inside `tmux -L hf` (the default tmux server still has the write stall; REVISION1_PLAN §8).
- `timings.txt` gives the run times for Table 1.

### C4.2 Runs

| # | Project | Command | Fills |
|---|---|---|---|
| 1 | `hf_1000g` | `run hf_1000g compare_gwases compare_gwases` | LDSC h² (liability) and rg for every pair; heterogeneity; expected vs observed (m1, R2-2, evaluation section) |
| 2 | `hf_1000g` | `run hf_1000g finemap finemap` | Fine-mapping HFall, 1000G LD, L auto (= 1) (M2, m11) |
| 3 | `hf_1000g` | `run hf_1000g coloc coloc` | Coloc HFall / ni-HFrEF / ni-HFpEF with N_eff; Fig. 2d (m5, m6) |
| 4–6 | `hf_1000g` | `for d in eqtlgen brainmeta ukb_pqtl; do run hf_1000g qtl_mr qtl_mr_$d; done` | MR-QTL (Fig. 2a) |
| 7 | `hf_1000g` | `run hf_1000g disease_progression disease_progression_slopehunter` | Disease progression, SlopeHunter (B27) |
| 8 | `hf_cwls` | `reuse hf_1000g hf_cwls; run hf_cwls disease_progression disease_progression_cwls` | Same with CWLS, for comparison |
| 9 | `hf_1000g_L10` | `reuse hf_1000g hf_1000g_L10; run hf_1000g_L10 finemap finemap_L10` | L = 10 vs L = 1 (M2, m11) |
| 10 | `hf_pan_ukb` | `run hf_pan_ukb finemap finemap_pan_ukb` | Fine-mapping with Pan-UKB LD (EUR N ≈ 420k) |
| 11 | `hf_pan_ukb` | `run hf_pan_ukb coloc coloc_pan_ukb` | Coloc with Pan-UKB LD |
| 12 | `hf_pan_ukb_L10` | `reuse hf_pan_ukb hf_pan_ukb_L10; run hf_pan_ukb_L10 finemap finemap_pan_ukb_L10` | L = 10 with Pan-UKB. N / N_REF ≈ 4, so `auto` already gives L = 10 here: run 10 and run 12 should agree, which is a useful check. |
| 13 | `hf_multi` | `run hf_multi finemap finemap_multi_ancestry` | MultiSuSiE HFall EUR + BBJ EAS, 1000G (m9) |
| 14 | `hf_multi_pan_ukb` | `run hf_multi_pan_ukb finemap finemap_multi_ancestry_pan_ukb` | MultiSuSiE with Pan-UKB (EAS N ≈ 2.7k) |

When everything has finished, restore the `.env`: `cp .env.backup .env`.

**Things to watch:**
- **Run 1 (LDSC):** `results/ldsc/results_EUR_h2.tsv` and `results_EUR_rg.tsv` must have three rows each.
  - The run mixes binary traits, so LDSC 1.0.1 gets `nan` prevalences only if a trait has no prevalence; here all three have one.
  - If `run_ldsc.sh` fails on the `--samp-prev`/`--pop-prev` lists, send me the log.
- **Strand check:** each standardised GWAS now has `<prefix>_std.allele_qc.tsv` (strand flips, palindromic variants, EAF discordance). Keep these for Supplementary Note X.
- **Build conversion QC:** `<prefix>_std.liftover_qc.tsv` gives the counts for the build-conversion sentence (R2-1 point 4).
- **Fine-mapping (runs 2, 9, 10, 12):** `results/finemap/FORMAT-METAL_Pheno1_EUR/` has `finemap_regions.tsv` (L, `LD_MISMATCH_S`, `N_LD_OUTLIERS`, `LD_FLAG`, MHC) and `finemap_credible_sets.tsv` (`LD_FLAG` per credible set). The HTML report is `results/finemap/result_finemap_FORMAT-METAL_Pheno1_EUR.html`.
- **CDKN1A (audit B13):** check the locus that had five credible sets in the submitted run. Compare its `LD_FLAG` and `N_CS` across runs 2, 9, 10 and 12.
- **Memory and time:** the LD-mismatch check does an eigen decomposition per region (once per ancestry in MultiSuSiE). Regions of 10,000+ variants take minutes each. If `run_finemapping` runs out of memory, lower its parallelism by setting `SLURM_MEM_PER_NODE` in the environment, or raise the rule's `mem`.
- **Disease progression:** check that the report says SlopeHunter for run 7 and CWLS for run 8.

## C5. LDSC benchmark (R1-M1)

**What it compares.** h² and rg for the same standardised GWAS (from run 1) with three LD score sets:
- Pan-UKB (the new default);
- the faulty v1.1.0 scores (`LDSCORE/b37_dbsnp156`, for the letter's "what went wrong");
- the Alkes Price group's 1000G Phase 3 scores (EUR, and EAS if wanted).

**Getting the Alkes Price scores.** They are in a requester-pays bucket, so a billing project is needed:

```bash
gsutil -u <billing-project> ls gs://broad-alkesgroup-public-requester-pays/LDSCORE/        # find the 1000G Phase 3 EUR weights / LD scores
gsutil -u <billing-project> cp gs://broad-alkesgroup-public-requester-pays/LDSCORE/<file>.tgz $WORK_DIR/
# unpack so that $WORK_DIR/ldscores_alkes/EUR/{1..22}.l2.ldscore.gz exists, with w_hm3.snplist in $WORK_DIR/ldscores_alkes/
```

**Running it.** `run_ldsc.sh` now:
- takes study types, prevalences and labels;
- runs every pair;
- writes `<prefix>_h2.tsv` and `<prefix>_rg.tsv`.

```bash
WORK_DIR=/local-scratch/data/genehackman_grch38; SIF=/home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif
STD=$OUT/hf_1000g/data/gwas
G=$STD/FORMAT-METAL_Pheno1_EUR_std.tsv.gz,$STD/FORMAT-METAL_Pheno3_EUR_std.tsv.gz,$STD/FORMAT-METAL_Pheno4_EUR_std.tsv.gz
NS=<Neff1>,<Neff3>,<Neff4>                 # fallback only: the standardised files have N_EFF, which is used
TYPES=binary,binary,binary; SAMP=NA,NA,NA; POP=<K1>,<K3>,<K4>; LABELS=HFall,ni-HFrEF,ni-HFpEF
mkdir -p $WORK_DIR/ldsc_benchmark
for set in pan_ukb:/home/wt23152/pipeline_data/LDSCORE/pan_ukb \
           v1.1.0:/home/wt23152/pipeline_data/LDSCORE/b37_dbsnp156 \
           alkes:$WORK_DIR/ldscores_alkes; do
  name=${set%%:*}; dir=${set#*:}
  apptainer exec -B /local-scratch,/home/wt23152 --env DATA_DIR=$WORK_DIR/ldsc_benchmark/$name $SIF \
    bash scripts/run_ldsc.sh "$G" "$NS" EUR "$WORK_DIR/ldsc_benchmark/$name/results_EUR.log" "$dir" \
    "$TYPES" "$SAMP" "$POP" "$LABELS"
done
head $WORK_DIR/ldsc_benchmark/*/results_EUR_{h2,rg}.tsv
```

Also look at the old v1.1.0 LD score build logs (`LDSCORE/b37_dbsnp156/<ANC>/*.log`) for the `--ld-wind-kb` flag; the letter asks for that confirmation.

**Done when** each set has an h² table (three traits) and an rg table (three pairs).

## C6. MungeSumstats concordance (R1-M3, R2-1)

The script is `scripts/validation/compare_with_mungesumstats.R`, and `scripts/validation/README.md` has the details.
- **What it does.** It runs both tools on the same file to GRCh38 with rsIDs, matches variants on position and allele pair, and reports the share of identical rsIDs, concordant effect directions and the correlation of effects.
- **`--chrpos_only`** rewrites the input so variants are identified only by `chr:pos:a1:a2`, which is the case Reviewer 1 said GeneHackman failed on and MungeSumstats handled.

**1. Install MungeSumstats into a user library for the container's R.** This is about 6 GB and needs internet.

```bash
SIF=/home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif
export R_LIBS_USER=/local-scratch/data/genehackman_validation/Rlib; mkdir -p $R_LIBS_USER
apptainer exec --env R_LIBS_USER=$R_LIBS_USER -B /local-scratch $SIF Rscript -e '
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP155.GRCh37", "SNPlocs.Hsapiens.dbSNP155.GRCh38",
                         "BSgenome.Hsapiens.1000genomes.hs37d5", "BSgenome.Hsapiens.NCBI.GRCh38"), update = FALSE, ask = FALSE)'
```

If Bioconductor refuses because of the R version in the container, install the Bioconductor release that matches it (`BiocManager::install(version = "<x.y>")`).

**2. Run the four comparisons** in `tmux -L mss`. Each takes a few hours and up to about 100 GB of RAM for a 10-million-variant file.

```bash
cd /home/wt23152/GeneHackman/scripts
V=/local-scratch/data/genehackman_validation; mkdir -p $V
mss() { apptainer exec --env R_LIBS_USER=$R_LIBS_USER,PIPELINE_DATA_DIR=/home/wt23152/pipeline_data \
          -B /local-scratch,/home/wt23152 $SIF Rscript validation/compare_with_mungesumstats.R --threads 16 "$@"; }
HERMES="SNP=rsID,CHR=chr,BP=pos_b37,EA=A1,OA=A2,EAF=A1_freq,BETA=A1_beta,SE=se,P=pval"
mss --gwas ../hf_data/FORMAT-METAL_Pheno1_EUR.tsv.gz --columns "$HERMES" --N <HFall N> \
    --input_build GRCh37 --output_build GRCh38 --output_prefix $V/hfall
mss --gwas ../hf_data/hf_japan_GCST90668009.tsv.gz --input_build <BBJ build> --output_build GRCh38 --output_prefix $V/bbj
# Reviewer 1's case: the same files identified only by chr:pos:a1:a2
mss --gwas ../hf_data/hf_japan_GCST90668009.tsv.gz --chrpos_only --input_build <BBJ build> --output_build GRCh38 \
    --output_prefix $V/bbj_chrpos
mss --gwas ../hf_data/FORMAT-METAL_Pheno1_EUR.tsv.gz --columns "$HERMES" --N <HFall N> --chrpos_only \
    --input_build GRCh37 --output_build GRCh38 --output_prefix $V/hfall_chrpos
```

**3. Check before reporting:**
- **Reviewer's case.** The `*_chrpos` runs should give (nearly) the same `variants_matched`, `rsid_identical_pct` and `effect_direction_concordant_pct` as the runs without it. If GeneHackman assigns no or few rsIDs there, the reviewer's bug is **not** fixed, and I need the log.
- **Discordant variants.** Read `<prefix>_discordant.tsv.gz`. Expected differences are:
  - multi-allelic sites (both tools pick one rsID);
  - palindromic variants (GeneHackman resolves them against the LD reference; MungeSumstats against the genome);
  - variants removed in build conversion by one tool and not the other.
  
  Anything else, for example many sign flips at non-palindromic SNVs, is a bug to look at before the numbers go in the letter.
- **Script failures.** The script itself hasn't been run end to end, because my sandbox couldn't reach the reference data. If MungeSumstats' column guessing trips over the rewritten input (`A1`/`A2`/`FRQ`), send me the error.

## C7. Release

These steps come from RUNBOOK §7 and §9:
1. **Upload 1.2.0.**
   ```bash
   gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r $WORK_DIR/bundle_1.2.0/ gs://genehackman/1.2.0/
   ```
   Before uploading the QTL bucket, decide on the eQTLGen, BrainMeta and UKB-PPP terms (REVISION1_PLAN §12 doesn't cover them).
2. **End-to-end tests** with the published data, installed the way users will install it (RUNBOOK §9), then `./tests/e2e_tests/run_test_pipelines.sh`. The finemap runs now also produce the HTML report, which is a quick check that `finemap.Rmd` renders in the container.
3. **Release.** Commit `tests/testing_complete.txt`, tag `1.2.0`, and record the image digest from C1.

## C8. What to send back, and what it fills

Send the files below. I'll turn them into the supplementary tables and fill the placeholders in `Manuscript - tracked 2026-11-06.docx` and `GeneHackman Response to Reviewers - tracked 2026-11-06.docx`.

| Placeholder (letter / manuscript) | File(s) |
|---|---|
| LD score benchmark, "differed by at most [X] SEs", Table S (M1) | `$WORK_DIR/ldsc_benchmark/*/results_EUR_{h2,rg}.tsv` |
| Liability-scale h², rg between subtypes, cross-trait intercepts (m1, m6, evaluation section) | `$OUT/hf_1000g/results/ldsc/results_EUR_{h2,rg}.tsv` |
| Fine-mapping: regions with LD mismatch, credible sets under L = 1 vs 10, 1000G vs Pan-UKB (M2, m11, evaluation section) | `finemap_regions.tsv` and `finemap_credible_sets.tsv` from runs 2, 9, 10, 12; the HTML reports |
| Colocalized signals with N_eff; Fig. 2d (m5) | `$OUT/hf_1000g/results/coloc/coloc_results.tsv` (with `ld_check`), the LocusZoom plots; run 11 for comparison |
| MR-QTL and disease progression sentences (evaluation section) | the `results/mr/` and `results/collider_bias/` outputs and reports of runs 4–8 |
| MungeSumstats concordance (M3, R2-1 point 3, manuscript testing sentence) | `$V/*_summary.tsv` (and a few lines of `*_discordant.tsv.gz` if anything looks odd) |
| Build conversion "agree with dbSNP for [X]% of [N]" (R2-1 point 4) | `validation_report.md` (REVISION1_PLAN §9.4) and the `*_std.liftover_qc.tsv` files |
| Strand / palindromic counts (Supplementary Note X) | `*_std.allele_qc.tsv` |
| Table 1 run times, Table S of run parameters | `$OUT/timings.txt`; the filled-in `paper/revision1_hf/*.yaml` |
| Dataset table (m9, m17): sources, N, N_eff, build, pipelines | the filled-in configs and C3.2 |
| Image tag/digest, reference data release | C1, C7 |

Decisions that are still yours (also flagged as comments in the documents):
- **m9:** BBJ in the GWAS comparison (option A) or a text fix (option B).
- **Disease progression:** report SlopeHunter, CWLS or both.
- **v1.1.0 LD scores:** keep `ldsc.ld_scores: 1000genomes_b37` (the faulty v1.1.0 set) selectable, or remove it.
- **QTL terms:** the terms of the QTL bucket's data.
