# GeneHackman revision 1, plan 2: run everything on ieup1

**6 November 2026.** Follows `REVISION1_PLAN.md` (reference data, parts A1–B2, all done) and `REVISION1_GAP_ANALYSIS.md` (what each reviewer point needs). This plan covers what still has to run on ieup1:

- the container with the new code;
- the 1.2.0 release;
- the heart-failure (HF) re-analysis that fills the yellow placeholders in the tracked manuscript and response letter;
- the LDSC benchmark;
- the MungeSumstats comparison, including Reviewer 1's `chr:pos:a1:a2` case.

Each part can be started on its own ("start part C4"), in this order:

| Part | What | Needs | Rough time | Status |
|---|---|---|---|---|
| **C0** | Get the code onto ieup1 | Branch pushed, or a bundle | 10 min | Done |
| **C1** | Container 1.2.0 (adds EnsDb v86): build, tag, push, check | C0 | 1–2 h (build) | Done (7 Nov) |
| **C2** | Unit tests, Python tests and dry runs on ieup1 | C0, C1 | 30 min | Done (7 Nov) |
| **C3** | HF data in place; fill in case and control counts | C0 | 30 min plus reading | Done (7–9 Nov) |
| **C4** | HF re-analysis: 14 pipeline runs | C1–C3 | 1–2 days of machine time | Done (9 Nov) |
| **C5** | LDSC benchmark (R1-M1) | C4 run 1 (standardised GWAS) | 1 h | Done (10 Nov) |
| **C6** | MungeSumstats concordance (R1-M3, R2-1), including Reviewer 1's case | C1, C3 | 1 day of machine time | Done (10 Nov) |
| **C7** | Release: upload 1.2.0, end-to-end tests, `testing_complete.txt` | C1, C2 | 1 day | Not started |
| **C8** | Send the results back for the documents | C4–C6 | – | Package ready (10 Nov) |

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

## C1. Container 1.2.0: build, tag and push

**Why.** The new code doesn't need a new image: the pipelines bind the repository's `R/` and `scripts/` into the container, so code changes take effect straight away. MungeSumstats isn't part of the pipeline and never goes in the image; it's installed separately, only for the validation script (C6). The image is rebuilt for two reasons:

- **The missing package.** GRCh38 gene tracks in LocusZoom (coloc, QTL MR) need `EnsDb.Hsapiens.v86`. It was added to `docker/requirements.R` on 1 Nov (`0e6b977`), and no image has been built since.
- **A citable tag.** The paper should cite 1.2.0 (audit B32). `DESCRIPTION` is now `Version: 1.2.0`, so `run_pipeline.sh` looks for `genehackman_1.2.0.sif` and pulls `docker://mrcieu/genehackman:1.2.0` if it isn't there. Push the image **before** running anything.

**`EnsDb.Hsapiens.v75` stays.** It's the GRCh37 gene annotation, used when a run sets `output.build: GRCh37`. That option is still supported (the letter says so, and the `*_grch37` end-to-end tests use it). It adds well under 100 MB to the image. Remove it only if GRCh37 output is dropped altogether.

1. **Build and push.** Run this from a machine with Docker and push access to `mrcieu`. `--platform linux/amd64` is needed on a Mac, and the build takes about 1–2 h there.
   ```bash
   cd GeneHackman && git checkout feature/1.1.1-reviewer-feedback
   docker build --platform linux/amd64 -f docker/Dockerfile -t mrcieu/genehackman:1.2.0 -t mrcieu/genehackman:develop .
   docker push mrcieu/genehackman:1.2.0 && docker push mrcieu/genehackman:develop
   docker inspect --format '{{index .RepoDigests 0}}' mrcieu/genehackman:1.2.0   # the digest to cite in the paper
   ```
2. **SIF on ieup1.** In `.env`, unset `DOCKER_VERSION` (or set it to `1.2.0`); `run_pipeline.sh` then builds the SIF on first use. To build it by hand:
   ```bash
   apptainer build /home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif docker://mrcieu/genehackman:1.2.0
   ```
3. **Check it.**
   ```bash
   SIF=/home/wt23152/pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif
   apptainer exec $SIF Rscript -e 'for (p in c("EnsDb.Hsapiens.v86","EnsDb.Hsapiens.v75","susieR","coloc")) cat(p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"), "\n"); cat("kriging_rss:", exists("kriging_rss", asNamespace("susieR")), "\n")'
   apptainer exec $SIF python -c "import MultiSuSiE, scipy; print('python ok')"
   ```
   Expect both EnsDb packages, susieR 0.12 or later with `kriging_rss: TRUE`, and `python ok`.

**C1 status (7 Nov 2026): done.**
- **Docker build.** The first build failed: `hgdownload.cse.ucsc.edu` serves a certificate for another hostname. The Dockerfile now downloads liftOver from `hgdownload.soe.ucsc.edu`.
- **Image.** `mrcieu/genehackman:1.2.0` was pushed 7 Nov, digest `sha256:12c536aad1e6830355ea70371ca1779b59559d6d989e5f99fc6d837e1713d88c`. This is the digest to cite.
- **SIF.** Built at `pipeline_data/genomic_data/pipeline/genehackman_1.2.0.sif`, and `.env` now has `DOCKER_VERSION=1.2.0`.
- **Package checks.** EnsDb.Hsapiens.v86 and v75 are both 2.99.0, susieR is 0.14.2 (`kriging_rss: TRUE`), coloc is 5.2.3, and MultiSuSiE and scipy import.

## C2. Tests and dry runs on ieup1

**C2 status (7 Nov 2026): done.** Logs are in `/local-scratch/data/genehackman_grch38/c2_logs_1.2.0/`.
- **R tests.** 145 test blocks, 0 failed, 0 errors, 1 skipped (the liftover chain test, which needs `GH_LIFTOVER_CHAIN_DIR`).
- **Python tests.** 20 tests, OK.
- **Dry runs.** All 14 HF configs plan without error (5–13 jobs each) against the 1.2.0 SIF, and the HERMES inputs map `N_CASES=N_case` and `N=N_total`.

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

The HERMES and BBJ files are already on ieup1 in `~/genehackman_project_dir` (`FORMAT-METAL_Pheno{1,3,4}_EUR.tsv.gz`, `hf_japan_GCST90668009.tsv.gz`). The configs in `paper/revision1_hf/` use the relative path `hf_data/…`, so the supplementary YAMLs have no local paths (Reviewer 1 m20). Link that folder into the repository root (`hf_data` is git-ignored):

```bash
cd /home/wt23152/GeneHackman && ln -sfn ~/genehackman_project_dir hf_data
ls hf_data/
```

**Missing:** `hf_progression_allcause_death.autosomes.gz`, the subsequent-event GWAS for disease progression (runs 7 and 8). Put it in `~/genehackman_project_dir` if it's elsewhere on ieup1; otherwise copy it from BluePebble (`/user/work/wt23152/020/`). Everything else can run without it.

`FORMAT-METAL_Pheno1_ALL` and `Pheno5-DCM` aren't used in the paper's analyses.

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

### C3 status (7 Nov 2026): done, except the prevalences

- **Builds.** BBJ and the progression GWAS are GRCh37. Their chr22 SNVs match the GRCh37 1000G panel (EAS 88%, EUR 92%) against 0.1% for GRCh38. The `build:` lines record this.
- **HERMES per-variant counts.** The files have `N_case` and `N_total` but no controls column. The configs map `N_CASES: N_case` and `N: N_total`.
  - GeneHackman now takes controls = `N` − `N_CASES` (`derive_n_controls()` in `R/gwas_formatting.R`, with tests).
  - So N_EFF is per variant. Per-variant N_EFF ranges: HFall 44k–487k, ni-HFrEF 2.0k–19.8k, ni-HFpEF 1.5k–14.3k. Every variant has valid counts.
  - The study counts below are used only where a variant has no counts.
- **Study counts.** From Henry et al. 2025 (Nat Genet 57:815, doi:10.1038/s41588-024-02064-3), Supplementary Table 2, European rows:
  - HFall: 139,533 cases and 1,568,809 controls.
  - ni-HFrEF: 4,975 cases and 935,830 controls.
  - ni-HFpEF: 3,590 cases and 752,599 controls. The 3,841 in the draft letter is the all-ancestry count.
- **Progression GWAS.** From Sunderland et al. 2026 (J Am Coll Cardiol 87:2981; the GWAS is in GWAS Catalog GCST90827692/3):
  - 50,636 European HF patients, of whom 11,454 died, which leaves 39,182 survivors. The file's maximum per-variant n is 50,434.
  - The GWAS is a Cox model, so N_EFF from deaths and survivors is an approximation.
  - The per-variant `n` is not mapped, as decided.
- **Prevalences (decided 9 Nov).** As Henry et al., the meta-analysis case fraction: HFall 0.08168, ni-HFrEF 0.005288, ni-HFpEF 0.004747. See the C4 status.

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

**C4 status (9 Nov 2026): all 14 runs done.** They ran one after another from `/local-scratch/data/genehackman_hf_runs/run_c4.sh` on 7 Nov, 15:23–20:53, about 5.5 h. Run times and peak memory are in `timings.txt`; the largest was 39 GB, for disease progression. Two runs failed and were rerun on 9 Nov (`rerun_c4_failed.sh`, logs `*.rerun.log`) after these fixes:
- **Run 1, LDSC.** munge_sumstats counted `N`, `N_EFF`, `N_CASES` and `N_CONTROLS` all as sample sizes ("Found 2 different N columns"). `run_ldsc.sh` now gives it only RSID, EA, OA, P, BETA, EAF and the chosen N.
  - Also fixed: with prevalences given, LDSC 1.0.1 prints only the liability-scale h². So `H2_Z` was empty, and every trait was flagged low-z and every rg as unreliable. `ldsc_results_tables()` now takes the z-score from the liability-scale estimate. The h² table needs regenerating, which happens with the prevalence rerun.
- **Run 13, MultiSuSiE with 1000G.** PLINK's trailing blank lines came back from `use_gwas_positions_in_clumped()` as a row of NAs, and that row was written into every `.clumped` file. MultiSuSiE then read it as a lead with no position.
  - Fixes: the R writer drops that row; the Python and R readers skip it; the Python reader restores the CHR and BP types.
  - Now 55 of 55 loci are fine-mapped (Pan-UKB run: 54 of 54).
  - The single-ancestry fine-mapping outputs had no NA regions, so they're unaffected.
- **Heterogeneity step.** `calculate_heterogeneity.R` passed `labels` to `create_dir_for_files()`. It now passes only the output paths. The outputs were written anyway.

**Prevalences (9 Nov 2026).** As Henry et al., the population prevalence is the meta-analysis case fraction: HFall 0.08168, ni-HFrEF 0.005288, ni-HFpEF 0.004747.
- Only LDSC (and the run 1 report) reran; dry runs of the other 13 runs showed nothing to do.
- Liability h²: HFall 3.93% (0.19), ni-HFrEF 7.80% (1.29), ni-HFpEF 4.18% (1.69, z = 2.5, flagged).

**C8 package (updated 10 Nov 2026 with C5 and C6).** Everything for the documents is in `/local-scratch/data/genehackman_hf_runs/for_documents/`, also as `genehackman_hf_results_2026-11-10.tar.gz` (which replaces the 9 Nov tarball).
- It contains the tables, reports, figures, QC files and configs.
- `README.md` gives the numbers by reviewer point, with caveats.
- `results_summary_tables.md` is generated by `summarise_results.py`.

**Earlier LDSC result, with the placeholder prevalences (0.02, 0.01, 0.01).**
- Liability h²: HFall 2.6% (SE 0.13%), ni-HFrEF 9.0% (1.5%), ni-HFpEF 4.9% (2.0%).
- rg: HFall–ni-HFrEF 0.67 (0.06), HFall–ni-HFpEF 0.80 (0.14), ni-HFrEF–ni-HFpEF 0.40 (0.18).
- Henry et al. (LDAK SumHer, not LDSC): h² 5.4%, 11.8% and 1.8%; rg ni-HFrEF–ni-HFpEF 0.42 (0.18).
- Henry et al. converted to the liability scale with the population prevalence set to the meta-analysis case fraction, and a sample prevalence of 0.5 with effective N. The EUR case fractions are 0.0817, 0.0053 and 0.0047; those prevalences would give liability h² of about 3.9%, 7.8% and 4.2%.

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

**C5 status (10 Nov 2026): done.** Script `$WORK_DIR/run_c5_ldsc_benchmark.sh`; results in `$WORK_DIR/ldsc_benchmark/<set>/results_EUR_{h2,rg}.tsv` and `comparison_with_alkes.tsv`.
- **Alkes Price set.** `eur_w_ld_chr` (1000G EUR, r² over HapMap3 SNPs, HLA removed; the LDSC tutorial's scores), from the requester-pays bucket with `-u brc-ieu`.
  - The plan named `1000G_Phase3_weights_hm3_no_MHC`, but that set has no `.l2.M` files, so it can only serve as weights.
- **Inputs.** The same standardised GWAS as C4 run 1, the same per-variant N_EFF, and the same prevalences (Henry case fractions).
- **Pan-UKB vs Alkes.**
  - Liability h² and rg differ by at most 0.75 Alkes SEs (HFall h² 3.93% vs 4.08%; rg 0.673 vs 0.684, 0.800 vs 0.774, 0.396 vs 0.340).
  - Including the intercepts, the maximum is 1.6 SEs (HFall intercept 1.038 vs 1.022).
- **v1.1.0 vs Alkes.**
  - HFall h² 11.9% vs 4.1% (39 Alkes SEs); ni-HFrEF h² 22.7% vs 8.5%.
  - HFall intercept 1.121 vs 1.022.
  - rg HFall–ni-HFpEF 1.22 (out of bounds).
- **`--ld-wind-kb 1.0`** is in all 110 v1.1.0 LD score logs (5 ancestries × 22 chromosomes, `LDSCORE/b37_dbsnp156/<ANC>/<chr>.log`). That is a 1 kb window instead of the usual 1 cM.
  - The short window makes the LD scores far too small, which inflates h² about threefold. rg is less affected, because both traits are biased the same way.
- **For the letter (M1).** "Differed by at most [X] SEs" is 0.75 for h² and rg, or 1.6 including intercepts.
  - The tables are in `for_documents/tables/ldsc_benchmark/`, and the write-up is in `for_documents/README.md` §6.

### Comparison with Henry et al. 2025 (m1), checked 10 Nov 2026

For the letter, present this positively: GeneHackman reproduces the published LDSC results, and the one clear difference comes from Henry et al. using a different heritability method.
- **Genetic correlations: like for like.** Henry et al. estimated rg with LDSC (Supplementary Table 19 says "estimated using LDSC"). All three subtype pairs agree to within 0.26 combined SEs:
  - HFall–ni-HFrEF 0.673 vs 0.658;
  - HFall–ni-HFpEF 0.800 vs 0.851;
  - ni-HFrEF–ni-HFpEF 0.396 vs 0.417.
  - That holds even with a different LD reference (Pan-UKB scores, against their 10,000 UKB Europeans) and per-variant N_eff.
- **LDSC intercept.** HF: Henry 1.01 (main text, no SE). GeneHackman 1.038 with Pan-UKB, 1.022 with the Alkes Price scores. Both mean polygenicity, not confounding.
- **h²: a different method.** Henry et al. used LDAK SumHer (BLD-LDAK, LDAK-Thin) with the same liability conversion.
  - HFall: 3.9% vs 5.4% (5.3 combined SEs; the only clear difference).
  - ni-HFrEF: 7.8% vs 11.8% (1.4 SE).
  - ni-HFpEF: 4.2% vs 1.8% (1.1 SE; both imprecise).
  - The conclusions are the same: modest common-variant heritability, with ni-HFrEF the most heritable.
  - The LDSC estimates are 66–73% of LDAK's for HFall and ni-HFrEF, which fits the two methods' different heritability models.
- **To do: find a general LDSC vs LDAK comparison to cite** for that last point, so the difference is explained by the literature rather than asserted. Candidates, not yet verified:
  - Speed et al. 2017 Nat Genet, "Reevaluation of SNP heritability in complex human traits";
  - Speed & Balding 2019 Nat Genet, "SumHer better estimates the SNP heritability of complex traits from summary statistics";
  - Speed et al. 2020 Nat Genet, "Evaluating and improving heritability models for complex traits".
  - Check what each actually shows before citing. If none fits, say only that the methods use different heritability models.
- **Written up.** In `for_documents/README.md` §3 (m1), with source values in `for_documents/tables/henry_2025_comparison.tsv`.

## C6. MungeSumstats concordance (R1-M3, R2-1)

**C6 status (10 Nov 2026): done.** Everything is in `/local-scratch/data/genehackman_validation/`. Results and logs per comparison:
- `<name>_summary.tsv`
- `<name>_discordant.tsv.gz`
- `<name>.log`
- `progress.txt` (times and peak memory)

The script is `run_c6_comparisons.sh`, run in `tmux -L mss`. Order: two chr22 tests on BBJ, then `hfall`, `bbj`, `bbj_chrpos`, `hfall_chrpos`, one after another, 02:44–04:48. Each full comparison took 27–33 min and used 50–65 GB of memory.

**Set-up.**
- **MungeSumstats 1.14.1** is installed in `Rlib/`. R only sees it with `--env R_LIBS_USER=$V/Rlib -B /local-scratch`.
- Bioconductor's download of `SNPlocs.Hsapiens.dbSNP155.GRCh38` broke off (HTTP/2 stream error). It was fetched with curl into `pkg/` and installed with `R CMD INSTALL`.
- **Fixed before the runs** (with tests):
  - `as_chrpos_identifiers()` returned a vector for data.table input, so the chr:pos runs couldn't have worked.
  - `change_column_names(remove_extra_columns = TRUE)` dropped every column when none was extra.
  - `compare_with_mungesumstats.R` now gives both tools only the mapped columns, so MungeSumstats can't misread HERMES's extra columns (`logP`, `N_case`, `p_het`).

**Results** (all GRCh37 → GRCh38, rsIDs assigned, variants matched on position and allele pair):

| Measure | HFall | HFall, chr:pos only | BBJ | BBJ, chr:pos only |
|---|---|---|---|---|
| Input variants | 9,871,737 | 9,871,737 | 7,968,587 | 7,968,587 |
| Kept by GeneHackman | 9,853,075 | 9,853,075 | 7,955,980 | 7,955,980 |
| Kept by MungeSumstats | 4,981,923 | 5,072,909 | 3,947,443 | 3,947,443 |
| Matched (kept by both) | 4,975,852 | 5,066,697 | 3,942,962 | 3,942,962 |
| Kept only by MungeSumstats | 6,071 | 6,212 | 4,481 | 4,481 |
| Identical rsID (of matched with both rsIDs) | 99.958% | 99.983% | 99.993% | 99.993% |
| rsID from MungeSumstats only | 57 | 137 | 89 | 89 |
| rsID from GeneHackman only | 0 | 0 | 0 | 0 |
| Concordant effect direction | 99.99996% | 99.99992% | 100% | 100% |
| Correlation of aligned effects | 0.999999 | 0.999997 | 1 | 1 |

**Reviewer 1's case is fixed.** With variants identified only by `chr:pos:a1:a2`, GeneHackman:
- keeps exactly the same variants as with full columns (9,853,075 and 7,955,980);
- assigns rsIDs to all matched variants;
- agrees with MungeSumstats as closely as in the normal runs.

For BBJ, which has no rsID column, the two modes give identical numbers. For HFall, only MungeSumstats changes between the modes: it validates the file's rsIDs when given them, and drops a different set.

**The discordant variants are the expected kinds.**
- **Different rsIDs at the same variant.** HFall 2,108 (42% palindromic); HFall chr:pos 843 (5% palindromic); BBJ 258 (15% palindromic).
  - dbSNP has several rsIDs at some positions (merged or allele-specific records), and the tools use different dbSNP builds (GeneHackman 156, MungeSumstats 155).
  - For HFall with rsIDs, MungeSumstats keeps the file's own rsID, while GeneHackman looks every variant up in dbSNP 156. That is why so many palindromic variants differ there and not in the chr:pos run.
- **Opposite effect signs.** 2 variants in HFall and 4 in HFall chr:pos, of about 5 million; none in BBJ.
  - Their effect sizes differ too (e.g. 0.0121 vs 0.0104 at 2:850299), so the tools took different input rows at a duplicated position. This is not an allele-alignment error.
  - One is a palindromic C/G variant (1:13091677), which the tools resolve differently.
- **Why MungeSumstats keeps about half the variants.**
  - With its defaults it removes variants at positions that are multi-allelic in dbSNP (`bi_allelic_filter`; on the chr22 test, 44,218 of 104,672).
  - It also removes indels that aren't in its SNP lookup.
  - GeneHackman keeps them (99.8% of the input survives build conversion). The concordance is computed on the variants both keep. The letter should say that MungeSumstats' defaults are stricter, not that GeneHackman retains spurious variants.

**For the letter and manuscript (M3, R2-1 point 3).** The table above, the "Reviewer 1's case is fixed" statement and the explanation of the discordant variants. All of it is also in `for_documents/README.md` §7.

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
3. **Release.** Commit `tests/testing_complete.txt`, tag the commit `1.2.0` (the image was tagged in C1), and put the image digest from C1 in the manuscript.

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
- **LDSC vs LDAK reference:** find and verify a general comparison to cite for the h² difference from Henry et al. (see C5, "Comparison with Henry et al.").
- **m9:** BBJ in the GWAS comparison (option A) or a text fix (option B).
- **Disease progression:** report SlopeHunter, CWLS or both.
- **v1.1.0 LD scores:** keep `ldsc.ld_scores: 1000genomes_b37` (the faulty v1.1.0 set) selectable, or remove it.
- **QTL terms:** the terms of the QTL bucket's data.
