# GeneHackman revision 1: plan for the next pieces of work (ieup1)

**4 November 2026; updated 5 November 2026, 15:30.** Follows on from `REVISION1_HANDOFF.md`. The work is in four parts that can be started one at a time ("start part A1"). Each part lists its own checks and when it is done.

| Part | What | Depends on | Status (5 Nov) |
|---|---|---|---|
| **A1** | Finish the Pan-UKB pack: allele-direction check, refresh GRCh38 keys, install as `genomic_data/pan_ukb_ld` | — | ✅ Done |
| **A2** | GRCh38 reference data: config fix, 1000G GRCh38 panel (PLINK 2 download), Pan-UKB LD scores, validation report | A1 for the validation report (A2.4) only | ✅ Data done and validated; code changes not yet committed |
| **B1** | QTL code: `eqtlgen`, `brainmeta`, `ukb_pqtl`, cis only; remove MetaBrain and trans | — | ✅ Done; committed as `4260371`, not pushed |
| **B2** | QTL data built from genotype-phenotype-map 1.0.0, GRCh38 then GRCh37 | B1 (layout), A1 (pack for liftover QC) | ✅ Done; `QTL_DATA_DIR` set in `.env` |

**What's left** (see §9 for the data and its sizes, §10 for the files to commit, §11 for the download script and §12 for the Pan-UKB redistribution decision):
1. Commit the download-script changes (§10), then push them along with `4260371` and `0e4b2e1`.
2. Upload 1.2.0 (runbook §7). The Pan-UKB redistribution question is decided (§12). The UKB-PPP, eQTLGen and BrainMeta terms for the QTL bucket haven't been checked.
3. The wrap-up in §6: the full end-to-end tests, and the `REVISION1_HANDOFF.md` updates.
4. Optional: free up space on `/local-scratch` (§9.3).

## 1. Context

The Pan-UKB LD pack has finished building. Andrew copied it to `/home/wt23152/pipeline_data/genomic_data/pan_ukb`. dbSNP 156 GRCh38 is at `…/genomic_data/dbsnp/b38_dbsnp156`.

Goal: finish handoff §3A and the data part of §3B, and replace the missing QTL data. When done, three things are true:

- `/home/wt23152/pipeline_data` (the `.env` `PIPELINE_DATA_DIR`) is a working **GRCh38 1.2.0 + GRCh37-extras** reference directory. Andrew chose the "live dir" option.
- GeneHackman's QTL data is rebuilt in GRCh38 from genotype-phenotype-map 1.0.0. That means eQTLGen, BrainMeta and UKB-PPP, all **cis only**. MetaBrain and trans are removed from GeneHackman.
- The unit tests and dry runs pass against the real data.

### What I verified on ieup1 (read-only)

**A1, build finished and validated: done.**
- Every ancestry has `full.done` and `validate.json` with `ok=True` and 0 problems.
- `MANIFEST.json` was written on 3 Nov.
- No build is running. So the runbook's worktree step isn't needed: work in `/home/wt23152/GeneHackman`, which is at `06fbf4b`, the same as origin.

**Pack location.**
- The copy is at `genomic_data/pan_ukb`, but the code reads `genomic_data/pan_ukb_ld`. See `R/constants.R:82`, `snakemake/util/constants.smk:120` and `50_assemble_bundle.sh`.
- The copy has a hand-added `README.md` that `$WORK_DIR/pack` doesn't have.

**A3, GRCh38 keys: still needed.**
- The annotation dates from 28 Oct, so it used the old liftover. EUR chr22 has no `round_trip` or `non_primary_contig` statuses.
- `hg38ToHg19.over.chain.gz` is missing from `$WORK_DIR/liftover_data/genomic_data/liftover/`. It is present in `pipeline_data/genomic_data/liftover/`.

**A2, allele direction:** not run. The GRCh37 1000G files and `.frq` are in `pipeline_data/genomic_data/1000genomes/b37_dbsnp156`.

**B1, dbSNP b38: OK.**
- The columns are the same as b37 (`RSID CHR BP REF ALT`), for chr1–24.
- rs4965031 sits 126 kb lower in b38 and on the inverted strand. The runbook says "~0.5 Mb", which is wrong.
- There is no `md5sums.txt`.

**Bug in `reference_data/grch38/config.sh:21`.** The `}` in `{CHR}` closes `${KG30X_FILE_PATTERN:-…}` early. The pattern becomes `…chr{CHR.filtered…vcf.gz}`, so every 30x VCF download would 404.

**1000G GRCh38 source.** Andrew asked for a prebuilt download. PLINK 2's resources page has `all_hg38` (pgen.zst 3.4 GB, plus pvar and psam), built from the **same NYGC 30x release**.
- It has 3,202 samples, a superpopulation column, and every one of the 2,504 phase-3 people.
- Multiallelic sites are split, and indels and SVs are included.
- No MAF filter has been applied.
- This replaces about 30 GB of per-chromosome VCFs.

**hg18ToHg38 chain:** missing from `pipeline_data/genomic_data/liftover/`. The UCSC URL returns 200.

**QTL data.**
- `gs://genehackman-qtl/1.1.0/` is empty. **Correction, 5 Nov:** it isn't. `gsutil du` shows 60.3 GB: `eqtlgen/` 29.3 GB and `metabrain/` 31.0 GB (see §9.2).
- genotype-phenotype-map (`/local-scratch/projects/genotype-phenotype-map/data/study/<study>/`) has `extracted_snps.tsv`, with lead chr, bp, log_p, ld_block and a `cis_trans` label.
- It also has `standardised/EUR_<chr>_<bp>.tsv.gz`, one LD-block region per lead. These are GRCh38, with columns `RSID CHR BP EA OA EAF BETA SE P SNP Z`. The keys are already GeneHackman's `CHR:BP_EA_OA`.
- N comes from `results/1.0.0/studies_processed.tsv.gz` (`sample_size`).
- Sources there are `eqtlgen` (eqtlgen-cis-*, 18,280), `brain_meta` (BrainMeta-cis-*, 16,729) and `ukb_ppp` (2,850 proteins, genome-wide, with cis and trans labelled per lead).

**Machine load.** genotype-phenotype-map R jobs are running. Use `THREADS=16`.

---

## 2. Part A1: finish the Pan-UKB pack

**Done when:** the A2 gate passes for all five ancestries; `update-keys` has run for all six; `/home/wt23152/pipeline_data/genomic_data/pan_ukb_ld` matches the new `MANIFEST.json`.

> **Status, 4 Nov 2026: done.**
> - **A1.1 allele-direction check (chr22 vs 1000G GRCh37) passed.** Results are in `/local-scratch/data/pan_ukb_ld/phase0/compare_<POP>_chr22.json`.
>
>   | Pop | EAF r | share closer if swapped | sign concordance, \|r\|>0.3 | r Pearson |
>   |---|---|---|---|---|
>   | EUR | 0.9987 | 0.014 | 0.99983 | 0.930 |
>   | EAS | 0.9976 | 0.016 | 0.99981 | 0.924 |
>   | CSA (vs SAS) | 0.9973 | 0.017 | 0.99985 | 0.916 |
>   | AFR | 0.9990 | 0.008 | 0.99948 | 0.875 |
>   | AMR | 0.9925 | 0.028 | 0.99747 | 0.852 |
>
> - **A1.2 `update-keys` ran for all six ancestries**, about 15 min each.
>   - Row counts are unchanged.
>   - `round_trip` is 0.0045–0.0059% of variants, lower than the runbook's guess of ~0.03%.
>   - `non_primary_contig` is 1,082–2,239 per ancestry.
>   - The new `MANIFEST.json` is at commit `06fbf4b`.
> - **A1.3 the pack is installed at `genomic_data/pan_ukb_ld`** (renamed from `pan_ukb`).
>   - Only the 132 `variants_chr*.tsv.gz` files and `MANIFEST.json` changed, and only they were copied.
>   - All 402 files match the new manifest's sha256s.
>   - `get_ld_matrix(..., ld_reference = "pan_ukb")` reads it correctly via `PIPELINE_DATA_DIR`.
> - **Also done:**
>   - `hg18ToHg38.over.chain.gz` is in `genomic_data/liftover/`.
>   - `dbsnp/b38_dbsnp156/md5sums.txt` is written.
>   - `hg38ToHg19.over.chain.gz` is copied into the build's `liftover_data`.

### A1.1 Handoff A2: allele-direction check (tmux, in parallel with A1.2)
- Run the handoff §3 A2 loop against the **copy**, `--pack /home/wt23152/pipeline_data/genomic_data/pan_ukb`.
  - A1.2 doesn't change its GRCh37 keys, EAF or LD.
  - Using the copy avoids racing `update-keys`.
- Settings:
  - `G1K=/home/wt23152/pipeline_data/genomic_data/1000genomes/b37_dbsnp156`;
  - a PLINK wrapper that binds `/local-scratch,/home/wt23152`;
  - output to `$WORK_DIR/phase0/compare_<POP>_chr22.json`.
- **Gate, for every ancestry:** EAF r ≈ 1, `share_closer_if_swapped` ≈ 0, and sign concordance ≥ 0.98 for |r| > 0.3. If any fails, **stop and report**. Everything downstream depends on the pack's orientation.

### A1.2 Handoff A3: refresh the pack's GRCh38 keys (tmux, ~1–2 h per ancestry)
```bash
cp /home/wt23152/pipeline_data/genomic_data/liftover/hg38ToHg19.over.chain.gz /local-scratch/data/pan_ukb_ld/liftover_data/genomic_data/liftover/
bash reference_data/pan_ukb_ld/run_build.sh update-keys      # all 6 ancestries, then finalise → new MANIFEST.json
```
In `logs/update_keys_<POP>.log`, check three things:
- `round_trip` ≈ 0.03%;
- `non_primary_contig` is present;
- the row counts are unchanged.

### A1.3 Install the pack in the live data dir (after A1.1 and A1.2)
```bash
cd /home/wt23152/pipeline_data/genomic_data && mv pan_ukb pan_ukb_ld
rsync -a /local-scratch/data/pan_ukb_ld/pack/ pan_ukb_ld/     # no --delete (keeps README.md)
```
- Check every file against the new `MANIFEST.json` sha256s with a one-off Python check.
- Add `hg18ToHg38.over.chain.gz` from UCSC to `genomic_data/liftover/`.
- Write `md5sums.txt` for dbsnp b38.

---

## 3. Part A2: GRCh38 reference data

**Done when:** `pipeline_data/genomic_data/1000genomes/b38_dbsnp156/` and `pipeline_data/LDSCORE/pan_ukb/` exist and pass their checks; `$WORK_DIR/validation_report.md` is written (needs A1 for the pack-key section).

> **Status, 5 Nov 2026 15:30: done.** All data is built, installed in `pipeline_data` and validated. Only the commit is left (§10).
> - **A2.1 code changes: made, not committed.** They are the `config.sh` brace fix, the `KG30X_SOURCE=plink2` path, the runbook, a single-core `chrpos_to_rsid()` path, the LD-score file choice, and (5 Nov) the `.frq` fix in `40_validate_grch38.R`. `local_config.sh` is written; its current contents are in A2.1 below.
> - **A2.4 validation report: done, all checks pass.** It is at `/local-scratch/data/genehackman_grch38/validation_report.md` (EUR and EAS, chr22). Per-population copies are `validation_report_{EUR,EAS}.md` and the logs are `validation_{EUR,EAS}.log`.
>
>   | Check | EUR | EAS | Expected |
>   |---|---|---|---|
>   | 1. Pack variants with a GRCh38 key | 129,840 / 130,488 | 109,563 / 110,201 | most |
>   | 1. Pack GRCh38 keys found in the 30x panel | 83.7% | 77.5% | most (UKB is imputed from other panels) |
>   | 1. Panel variants with MAF ≥ 5% found in the pack | 89.0% | 85.3% | — |
>   | 1. EAF r, pack vs panel (mean \|diff\|; n with \|diff\| > 0.2) | **0.9980** (0.013; 113) | **0.9961** (0.017; 130) | > 0.95, very few > 0.2 |
>   | 2. GRCh37 panel lifted with GeneHackman's liftover | 119,883 / 121,435 | 104,918 / 106,428 | — |
>   | 2. Same position as the 30x panel, shared rsIDs | 99.996% (4 of 113,570 differ) | 99.997% (3 of 98,725) | ~100% |
>   | 3. Panel continuity: r correlation, GRCh37 vs GRCh38 panel pairs | 0.9925 (1.05 M pairs) | 0.9911 (1.04 M pairs) | > 0.95 |
>   | 3. Sign agreement for \|r\| > 0.3 | 99.84% | 99.75% | ~100% |
>   | 4. Expected files present | 25 / 25 | 25 / 25 | all |
>
>   - The positions that differ in section 2 are 3–4 rsIDs, the same ones in both populations, moved by 20 bp to 1.1 kb. Most likely dbSNP re-mapped them.
>   - **Run against the live dir.** I set `BUNDLE_DIR=/home/wt23152/pipeline_data`, because `bundle_1.2.0` only gets the chains, dbSNP and the pack at runbook step 7. So section 4 checks the live dir. I also set `PLINK19=plink1.9`, because `local_config.sh`'s `PLINK19` is a wrapper that calls `apptainer` and fails inside the container. Both points are now noted in RUNBOOK §6.
>   - **Bug fixed in `40_validate_grch38.R`.** On the first run, section 1 gave EAF r = 0.003. `read_panel()` treated the `.frq` `MAF` as the frequency of the `.bim` A1, but PLINK 1.9 `--freq` reports the minor allele as A1. That allele differs from the `.bim` A1 for about 49% of rows, in both the GRCh37 and GRCh38 panels. The script now flips MAF where the two A1s differ. GeneHackman's own `.frq` reader (`R/data_conversions.R`, EAF fill) and the A1.1 check (`61_compare_1000g.py`) already match on the allele, so neither is affected.
> - **A2.3 LD scores: done.**
>   - **Source files:** the archive (`UKBB.ALL.ldscore/`) has, per ancestry, `UKBB.<POP>.rsid.l2.ldscore.gz` (SNP = rsID) and `UKBB.<POP>.l2.M`/`.M_5_50`, plus stratified 8LDMS and 25LDMS sets for S-LDSC. The importer now uses the `.rsid` single-annotation files.
>   - **Counts** (variants per ancestry, % in `w_hm3.snplist`): EUR 1,094,429 (98.4%), AFR 1,190,493 (90.7%), AMR 1,143,767 (95.4%), CSA 1,116,449 (97.1%), EAS 997,351 (97.0%), MID 1,184,445 (92.6%).
>   - **M:** M = 23.96 M and M_5_50 = 7.01 M for EUR. These are Pan-UKB's own counts; LDSC uses M_5_50 by default.
>   - **Installed** in `pipeline_data/LDSCORE/pan_ukb/` (md5-checked). `run_ldsc.sh` already reads this layout and maps SAS to CSA.
> - **A2.2 panel: done (5 Nov 14:42), installed** in `pipeline_data/genomic_data/1000genomes/b38_dbsnp156/` (sizes and md5 checked).
>   - **Source:** PLINK 2's `all_hg38` fileset (3,202 people, of whom all 2,504 phase 3 people are present).
>   - **Variants / people:** EUR 9,006,714 / 503; EAS 7,819,052 / 504; AFR 15,218,110 / 661; AMR 10,079,712 / 347; SAS 9,265,049 / 489. That's 3–5% more than each GRCh37 panel. 96.5% of the GRCh37 EUR rsIDs are present.
>   - **rsIDs:** 20,301,229 of 20,303,769 variants have a dbSNP 156 rsID; 72,180 rsIDs that name several variants were dropped.
>   - **EUR checks:** no duplicate or non-rs IDs, A1 always alphabetical, nothing below MAF 1%.
>   - **First attempt lost:** the overnight run kept its temp dir in `/dev/shm`. Its files were deleted at 06:15, probably by logind's `RemoveIPC`, and the merge then failed. The rerun used `$WORK_DIR/run_panel_a2.sh` (parallel ranged download, then the build), with the temp dir on `/local-scratch`.
> - **ieup1 write stall: the cause is cgroup-specific.**
>   - The ~122 GB of dirty pages that don't drain are all charged to `session-9570.scope`. That's the cgroup of the default tmux server (sessions `a1`, `b`, `full_ld_blocks`, `redo_svgs`).
>   - Every write from a process started in that tmux server is throttled in `balance_dirty_pages`. Writes from other sessions run at full speed.
>   - **Workaround:** start long jobs on a separate tmux server, `tmux -L <name>`, launched from a fresh login. A2.2 ran on `tmux -L a2`.
>   - `$WORK_DIR/1000g_30x` is now a real directory, and the `/dev/shm` copies are gone.

### A2.1 Small fixes (commit)
- **`reference_data/grch38/config.sh`**
  - Fix the brace bug: put the default in a plain variable first.
  - Add `KG30X_SOURCE` (`plink2` default, or `vcf`).
  - Add overridable `PLINK2_HG38_PGEN_URL`, `PLINK2_HG38_PVAR_URL` and `PLINK2_HG38_PSAM_URL` (the Dropbox links from the PLINK 2 resources page, 2022-08-04 Byrska-Bishop release).
- **`reference_data/grch38/10_build_1000g_panel.sh`**, `plink2` source:
  - Download the three files once into `$TMP/plink2_hg38/`.
  - Run `plink2 --zst-decompress` on the pgen. Name the pvar `all_hg38.pvar.zst` and the psam `all_hg38.psam`.
  - Per chromosome, use `--pfile $TMP/plink2_hg38/all_hg38 vzs` in place of `--vcf`. Every other flag and step stays as it is: `--keep` the 2,504, `--snps-only just-acgt --max-alleles 2`, IDs `@:#:$r:$a`, `--rm-dup`, MAF per superpopulation, `11_assign_rsids.R` from our dbSNP .fst, and A1 alphabetical.
  - Keep all 2,504 people, as the GRCh37 panel does. Don't apply the deg2 relatives list.
  - Split multiallelic records that share an rsID are already dropped by `11_assign_rsids.R`. Confirm this in the log.
- **`reference_data/grch38/RUNBOOK.md`**
  - Step 0: no worktree is needed once the build is done.
  - Step 1: fix the rs4965031 expectation.
  - Step 2: describe the plink2 source.
  - Step 4: point to the new QTL route (Parts B1 and B2).
  - Step 5: copying the chain from the old data dir is enough.
- **Untracked `reference_data/grch38/local_config.sh`.** The root `.gitignore` already covers `local_config.sh`. Its contents as of 5 Nov:
  ```bash
  export WORK_DIR=/local-scratch/data/genehackman_grch38
  export BUNDLE_DIR=$WORK_DIR/bundle_1.2.0
  export OLD_PIPELINE_DATA_DIR=/home/wt23152/pipeline_data
  export DBSNP_B38_DIR=/home/wt23152/pipeline_data/genomic_data/dbsnp/b38_dbsnp156
  export LD_PACK_DIR=/home/wt23152/pipeline_data/genomic_data/pan_ukb_ld
  export NEW_QTL_DATA_DIR=/local-scratch/data/genehackman_grch38/qtl_1.2.0
  export CHAIN_DIR=/home/wt23152/pipeline_data/genomic_data/liftover
  export R_CMD="apptainer exec --bind /local-scratch,/home/wt23152 /local-scratch/data/pan_ukb_ld/genehackman_develop.sif Rscript"
  export PLINK2=/usr/local/bin/plink2
  export PLINK19=/local-scratch/data/pan_ukb_ld/plink1.9   # host wrapper; use PLINK19=plink1.9 for scripts run via $R_CMD
  export PATH=/local-scratch/data/genehackman_grch38/env/bin:$PATH
  export THREADS=16
  ```
- **Tools env:** `conda env create -p /local-scratch/data/genehackman_grch38/env -f reference_data/grch38/environment.yml`.

### A2.2 Runbook step 2: GRCh38 1000G panel (tmux)
- Run `bash reference_data/grch38/10_build_1000g_panel.sh | tee $WORK_DIR/panel.log`. It writes to `$BUNDLE_DIR/genomic_data/1000genomes/b38_dbsnp156`.
- Check the counts:
  - people: EUR 503, EAS 504, AFR 661, AMR 347, SAS 489;
  - variants: EUR ≈ 8.5 M.
- Then `cp -a` it to `pipeline_data/genomic_data/1000genomes/b38_dbsnp156/`.

### A2.3 Runbook step 3: Pan-UKB LD scores
- Run `python reference_data/grch38/20_import_pan_ukb_ldscores.py | tee $WORK_DIR/ldscores.log`. The archive is 2.4 GB.
- Read the log for:
  - the archive layout;
  - whether SNPs are rsIDs or had to be mapped through the pack;
  - the `w_hm3` overlap.
- Then `cp -a $BUNDLE_DIR/LDSCORE/pan_ukb /home/wt23152/pipeline_data/LDSCORE/`.

### A2.4 Runbook step 6: validation report
- Run `$R_CMD reference_data/grch38/40_validate_grch38.R --pop EUR --chrom 22`, then EAS. This writes `$WORK_DIR/validation_report.md`.
- The expected values are in runbook §6. The pack's GRCh38 keys should match the native panel with EAF r > 0.95.

---

## 4. Part B1: QTL code (cis only, three datasets)

**Done when:** testthat, Python tests and the QTL dry runs pass; MetaBrain and trans are gone from code, docs and test inputs.

> **Status, 5 Nov 2026: done and committed** as `4260371` ("QTL MR: cis-only eQTLGen, BrainMeta and UKB-PPP from the genotype-phenotype map") on `feature/1.1.1-reviewer-feedback`. It is **not pushed**: the branch is 1 ahead of origin. The worktree `/home/wt23152/GeneHackman_qtl` and its branch `feature/1.2.0-qtl-cis` are gone.
> - **Tests:**
>   - testthat: 369 passed, 0 failed.
>   - Python unittest, in the develop container: 9 passed.
>   - Snakemake dry runs plan all three datasets, and reject `subcategory: trans` and `dataset: metabrain`.
> - **Smoke runs on the real GRCh38 QTL data:** `perform_mr_on_qtl_dataset()` on `test_data_small_hg38.tsv.gz`:
>   - tests 145 eQTLGen, 145 BrainMeta and 10 UKB-PPP exposures, one row each, with F-statistics;
>   - `load_qtl_gwas_for_mr_result()` loads the coloc region (3,000–5,000 SNPs) for the top result of each.
> - **Changes from the plan:**
>   - No per-dataset default N. Every `gwas/cis` file carries N, so `qtl.N` stays only as the fallback.
>   - `--ancestry` and `--subcategory` are removed from `run_mr_against_qtl_data.R`, not ignored.
>   - Docs updated: `snakemake/qtl_mr.md` and `README.md`, plus `reference_data/grch38/{README,RUNBOOK}.md` and `config.sh` for the new QTL route.

### B1.1 GeneHackman code: three cis-only datasets with one generic layout (commit)

**New layout.**
- `QTL_DATA_DIR/<dataset>/<build>/top_hits/cis.tsv.gz`, with the columns `SNP CHR BP EA OA EAF BETA SE P N RSID ENSEMBL_ID EXPOSURE`.
- `QTL_DATA_DIR/<dataset>/<build>/gwas/cis/<EXPOSURE>.tsv.gz`.
- This is the 1.1.0 eQTLGen convention, now used for every dataset.

**Changes by file.**
- **`R/constants.R`**
  - `qtl_datasets <- list(eqtlgen="eqtlgen", brainmeta="brainmeta", ukb_pqtl="ukb_pqtl")`.
  - Remove `pqtl_top_hits_dir`.
  - Add a per-dataset default N, used only as a fallback because the files carry N.
- **`R/mr.R`**
  - Replace `perform_mr_on_metabrain_datasets` and `perform_mr_on_eqtlgen_datasets` with a single `perform_mr_on_qtl_dataset(gwas_filename, dataset, exposures, results_output, build)`.
  - Keep `run_mr_on_qtl_data` as it is.
  - Also fix the handoff's `subset(mr_results, exposure %in% …)` bug when `exposures` is set, if it is in this function.
- **`R/coloc.R`**: make `load_qtl_gwas_for_mr_result()` generic (`gwas/<outcome>/<EXPOSURE>.tsv.gz`), and drop the metabrain branch.
- **`scripts/run_mr_against_qtl_data.R` and `scripts/coloc_of_mr_results.R`**: one generic call; drop `--ancestry` and `--subcategory`, or ignore them.
- **`snakemake/qtl_mr.smk`**
  - Validate `qtl.dataset` against the three datasets.
  - `qtl.subcategory` becomes optional, defaults to `cis`, and anything else errors with "only cis QTLs are supported". `qtl_name` stays `<dataset>_cis`.
  - Use the dataset's default N when `qtl.N` isn't given.
- **Docs and templates**
  - `snakemake/qtl_mr.md`: dataset table with sources, N and citations: eQTLGen (Võsa 2021), BrainMeta (Qi 2022), UKB-PPP (Sun 2023).
  - `snakemake/input_templates/qtl_mr.yaml` and the README QTL section.
- **Tests**
  - Update `tests/testthat/test_mr.R` and `test_reference_builds.R`.
  - Replace `qtl_mr_metabrain.yaml` with `qtl_mr_brainmeta.yaml`, and add `qtl_mr_ukb_pqtl.yaml`.
  - `tests/e2e_tests/run_test_pipelines.sh` runs all three when `QTL_DATA_DIR` is set.

---

## 5. Part B2: QTL data from genotype-phenotype-map

**Done when:** `qtl_1.2.0/<dataset>/{GRCh38,GRCh37}/` exist for all three datasets, the sanity checks below pass, and `.env` has `QTL_DATA_DIR` set.

> **Status, 5 Nov 2026: done.** The data is built and checked, and `.env` has `QTL_DATA_DIR=/local-scratch/data/genehackman_grch38/qtl_1.2.0`. In total it takes 9.0 GB (§9).
>
> **GRCh37 copy:** `30_lift_qtl_data.R` dropped 0.09% of variants for eQTLGen, 0.06% for BrainMeta and under 0.01% for UKB-PPP. Sizes are the same as GRCh38.
>
> **Sanity checks vs the old fixtures**, on GRCh37 positions:
> - **eQTLGen** (`some_eqtlgen_hits`, 394 genes):
>   - 338 genes have an instrument in the new data;
>   - 116 have the same lead SNP;
>   - the fixture lead is in the new coloc region for 237. For those, the BETA sign agrees 237/237, and BETA and −log10 P correlate at 1.000.
> - **UKB-PPP** (`ukb_pqtl_cis_only`, 1,162 genes):
>   - 953 genes have an instrument;
>   - 527 have the same lead;
>   - the fixture lead is in the region for 897. The sign agrees 897/897, with r = 1.000.
> - **Why leads differ:**
>   - ties at P = 0: 307 UKB-PPP genes, 60 eQTLGen;
>   - fixture leads outside the cis LD block: 94 eQTLGen. The old eQTLGen fixture is genome-wide, and its median distance to the new lead is 0.3 Mb, upper quartile 30 Mb;
>   - fixture leads that are not the strongest SNP in the block.
>   - In only one case is the new lead weaker than the fixture's lead, at P ≈ 1e-300.
>
> GRCh38 build details (log `qtl_from_gpm.log`; output in `/local-scratch/data/genehackman_grch38/qtl_1.2.0/<dataset>/GRCh38/`):
>
> | Dataset | Studies | With a cis lead | Exposures with an instrument (P < 5e-8) | Instruments with an rsID | gwas files | Size |
> |---|---|---|---|---|---|---|
> | eqtlgen | 18,280 | 16,917 | 14,630 | 14,625 | 16,904 | 1.7 GB |
> | brainmeta | 16,729 | 16,711 | 15,287 | 15,233 | 16,590 | 2.2 GB |
> | ukb_pqtl | 2,850 | 2,364 | 1,658 | 1,647 | 2,362 | 0.4 GB |
>
> - **One cis lead per study.** Every study has at most one: genotype-phenotype-map labels a single LD block per gene as cis. So every instrument is a single SNP and every MR test is a Wald ratio.
> - **eQTLGen:** 1,363 studies have no `extracted_snps.tsv`. Presumably no region reached genotype-phenotype-map's extraction threshold; not checked.
> - **UKB-PPP:** 323 cis leads were dropped because they are 1.0–4.8 Mb from the gene (median 1.6 Mb). All are on the gene's own chromosome. genotype-phenotype-map's cis label is wider than the usual ±1 Mb (Sun et al. 2023).
> - **No region file:** 13 eQTLGen, 121 BrainMeta and 2 UKB-PPP studies have a cis lead but no standardised region file. These are dropped.

### B2.1 Conversion script: `reference_data/grch38/31_qtl_from_gpm.R` (run via `$R_CMD`, multi-core)

**Inputs:** `GPM_STUDIES=/local-scratch/projects/genotype-phenotype-map/results/1.0.0/studies_processed.tsv.gz`. For each dataset (`eqtlgen` ← source `eqtlgen` with `study_name` starting `eqtlgen-cis-`; `brainmeta` ← `brain_meta`; `ukb_pqtl` ← `ukb_ppp`), per study:

1. **Cis leads.** Take the rows of `extracted_snps.tsv` with `cis_trans == "cis"`.
   - For UKB-PPP, also require the lead within ±1 Mb of the gene. The gene coordinates come from genotype-phenotype-map's `data/gene_ranges_hg38.tsv`. This double-checks its label.
2. **Coloc region.** Read the `standardised/EUR_<chr>_<bp>.tsv.gz` file of each cis lead.
   - Join them, removing duplicate SNPs, and add `N` = `sample_size`.
   - Write `gwas/cis/<EXPOSURE>.tsv.gz`.
3. **Top hits.**
   - For each cis lead with P < 5e-8, take the lead SNP's row from its region file.
   - Each lead is in a distinct LD block, so they are roughly independent, and IVW applies when a gene has more than one.
   - Fill `RSID` from dbSNP b38 with GeneHackman's own rsID lookup.
   - Append to `top_hits/cis.tsv.gz`.

**Exposure names and P = 0.**
- `EXPOSURE` is the gene symbol when it is unique within the dataset. Otherwise it is `<symbol>_<ENSG or OID>`, so file names are unique.
- Keep P = 0 rows (UKB-PPP `log_p = Inf`). MR and coloc use BETA and SE.

**Run.**
- Resumable: skip exposures already written.
- Write `$NEW_QTL_DATA_DIR/<dataset>/GRCh38/summary.tsv` with counts per dataset: studies, cis leads, top hits, and exposures with no cis lead.
- Expected size is roughly a few GB to ~15 GB in total.

### B2.2 GRCh37 copies (for `output.build: GRCh37`)
- Generalise `reference_data/grch38/30_lift_qtl_data.R` to take `--from` and `--to` builds. It currently lifts 37→38 only.
- Lift `<dataset>/GRCh38` to `<dataset>/GRCh37` with `convert_reference_build_via_liftover()`, using the same QC as a GWAS.
- `liftover_summary.tsv` should show a drop of 5% or less.
- Then set `QTL_DATA_DIR=/local-scratch/data/genehackman_grch38/qtl_1.2.0` in `.env`.

---

---

## 6. Validation, per part
- **A1:** the gates in A1.1 and A1.2, plus the MANIFEST checksum check after the rsync in A1.3.
- **A2:** the checks in A2.2, A2.3 and A2.4.
- **QTL sanity checks**
  - eQTLGen top hits vs the 1.1.0 fixture `tests/testthat/data/some_eqtlgen_hits.tsv.gz`: match genes and keys after lifting the fixture to 38. Expect the same lead or high LD, and the same BETA sign.
  - UKB-PPP vs `tests/testthat/data/ukb_pqtl_cis_only.tsv`: the same check.
  - Spot-check that a known strong cis-eQTL is present in all three datasets.
- **Tests**, all in the develop container:
  - testthat with `GH_LIFTOVER_CHAIN_DIR=/home/wt23152/pipeline_data/genomic_data/liftover`;
  - Python unittest;
  - Snakemake `-n` dry runs of every e2e input, including the three QTL ones, against the live `pipeline_data`;
  - smoke runs of `scripts/run_mr_against_qtl_data.R` and `coloc_of_mr_results.R` on `test_data_small.tsv.gz` for each dataset.
- **Wrap-up:** update `REVISION1_HANDOFF.md` (§1 data table, §3 status). Andrew needs the A1.1 JSONs, `validation_report.md`, the A2.3 log summary and the QTL `summary.tsv`.

## 7. Out of scope (later)
- Container rebuild with EnsDb v86, full e2e runs, and the 1.2.0 release (§3 D).
- Bundle assembly and upload to `gs://genehackman` and `gs://genehackman-qtl/1.2.0/`.
  - Pan-UKB pack and LD scores: decided, they can be re-hosted (§12).
  - Not yet checked: the terms for UKB-PPP (UKB terms), BrainMeta and eQTLGen summary statistics in the QTL bucket.
- LDSC benchmark (B8): needs the paths of the HF standardised GWAS.
- Decisions C1–C3, other §3 E work, §3 F, and the optional A4/A5.

## 8. Risks and open points
- ~~**A2 fails:** the pack's orientation is wrong.~~ Resolved: the A1.1 check and the A2.4 report both show the pack correctly oriented (EAF r ≥ 0.996).
- **The ieup1 write stall is still there.** The default tmux server's cgroup is still holding about 122 GB of dirty pages (§3, A2 status). Run long jobs under `tmux -L <name>` until that tmux server is restarted.
- **genotype-phenotype-map regions:**
  - **Coverage.** Its regions were only extracted for leads with p < 1.5e-4 and are LD-block-sized, so the coloc "gwas" files cover the cis lead's block, not the full ±1 Mb cis window. This is fine for coloc around the MR hit, but worth a methods sentence.
  - **Imputation.** These are pre-imputation, as Andrew intended.
- **Dropbox links** for the PLINK 2 files could move. They are overridable in `config.sh`, and `KG30X_SOURCE=vcf` stays as a fallback, now with the brace bug fixed.
- **Manuscript and response letter:** dropping MetaBrain changes the QTL datasets they describe.

## 9. Where the data is and how much space it takes (5 Nov 2026)

Sizes are in GB (10⁹ bytes): local ones from `du`, bucket ones from `gsutil du -s`.

### 9.1 Local data GeneHackman reads

**`/home/wt23152/pipeline_data`** (`.env` `PIPELINE_DATA_DIR`): **123.3 GB** in total.

| Path | GB | Build | Bucket it goes to |
|---|---|---|---|
| `genomic_data/pan_ukb_ld/` (6 ancestries, `MANIFEST.json`, `README.md`) | 74.50 | 37 + 38 keys | `genehackman/1.2.0` |
| `genomic_data/dbsnp/b38_dbsnp156/` | 11.66 | 38 | `genehackman/1.2.0` |
| `genomic_data/1000genomes/b38_dbsnp156/` (new, A2.2) | 10.92 | 38 | `genehackman/1.2.0` |
| `LDSCORE/pan_ukb/` (6 ancestries + `w_hm3.snplist`, A2.3) | 0.09 | rsID | `genehackman/1.2.0` |
| `genomic_data/liftover/` (5 chains + `liftOver`) | 0.01 | — | both |
| `genomic_data/gene_name_map.tsv` | 0.003 | — | `genehackman/1.2.0` |
| `genomic_data/dbsnp/b37_dbsnp156/` | 11.44 | 37 | stays in `genehackman/1.1.0` (`--grch37`) |
| `genomic_data/1000genomes/b37_dbsnp156/` | 10.54 | 37 | stays in `genehackman/1.1.0` (`--grch37`) |
| `LDSCORE/b37_dbsnp156/` | 1.09 | 37 | stays in `genehackman/1.1.0` (`--grch37`) |
| `genomic_data/pipeline/genehackman_1.1.0.sif` | 3.01 | — | not uploaded (1.1.0 container) |

**`/local-scratch/data/genehackman_grch38/qtl_1.2.0`** (`.env` `QTL_DATA_DIR`): **9.02 GB** in total.

| Dataset | GRCh38 (published) | GRCh37 (local only, not published) |
|---|---|---|
| eqtlgen | 1.81 | 1.81 |
| brainmeta | 2.30 | 2.30 |
| ukb_pqtl | 0.40 | 0.40 |
| **Total** | **4.51** | **4.51** |

### 9.2 The two buckets: now and after the 1.2.0 upload

Nothing has been uploaded yet. The Pan-UKB redistribution question is decided (§12).

| Bucket | Now | Adds for 1.2.0 | After (if 1.1.0 is kept) |
|---|---|---|---|
| `gs://genehackman` | **23.1 GB**: `1.1.0/` = 1000genomes b37 10.54, dbsnp b37 11.44, liftover 0.01, LDSCORE b37 1.08 | `1.2.0/` 97.2 GB (pack 74.5, dbSNP b38 11.7, panel b38 10.9, LD scores 0.09, chains, gene map) | **120.3 GB** |
| `gs://genehackman-qtl` | **60.3 GB**: `1.1.0/eqtlgen` 29.3, `1.1.0/metabrain` 31.0 | `1.2.0/` 4.4 GB (three datasets, GRCh38 only) | **64.7 GB** |

- **Decided 5 Nov: no GRCh37 copy in 1.2.0.** `1.2.0` is GRCh38 only. The GRCh37 files for `output.build: GRCh37` stay in `1.1.0/` and are fetched from there by `download_reference_data.sh --grch37`. The md5s GCS stores for the 1.1.0 objects match an independent decode for all 944 objects, and match the local files I spot-checked. The QTL data are GRCh38 only, so `qtl_mr` requires `output.build: GRCh38`.
- **What a user stores.** Users download only what they need (§11), not the whole bucket: a default EUR install is 19.7 GB (7.8 GB with `--no-dbsnp`), against 97 GB for all of 1.2.0.
- The QTL bucket's 1.1.0 data (eQTLGen genome-wide plus MetaBrain, 60 GB) is about 6.7 times the new cis-only 1.2.0 set.

### 9.3 Working space on `/local-scratch` (nothing deleted)

| Path | GB | What it is | Can it go? |
|---|---|---|---|
| `genehackman_grch38/` (all) | 47.8 | | |
| ↳ `qtl_1.2.0/` | 9.0 | live QTL data | **keep** |
| ↳ `1000g_30x/` | 27.0 | A2.2 temp files (`all_hg38` pgen, per-chromosome pgen/bed) | yes, once the panel has been in use for a while; rebuilding takes about 30 min |
| ↳ `bundle_1.2.0/` | 11.0 | panel + LD scores, already copied to `pipeline_data` | yes; step 7 refills it anyway |
| ↳ `env/`, `qtl_work/`, `1000g_30x_partial_download/` | 0.8 | tools env, B2 temp files, a stale partial download | `env` keep; others yes |
| `pan_ukb_ld/` (all) | 180.6 | | |
| ↳ `pack/` | 74.5 | the same pack as `pipeline_data/genomic_data/pan_ukb_ld` | yes, once you're sure the build won't be rerun (`LD_PACK_DIR` now points to `pipeline_data`) |
| ↳ `work/` | 72.8 | build intermediates | probably yes |
| ↳ `logs/` | 14.7 | build logs | could be compressed |
| ↳ `source/`, `raw/` | 13.1 | downloaded Pan-UKB source | keep while a rebuild is possible |
| ↳ `genehackman_develop.sif` | 3.0 | the container `R_CMD` uses | **keep** |
| ↳ `env/`, `phase0/` | 2.5 | tools env; A1.1 JSONs | **keep** |

Deleting `1000g_30x`, `bundle_1.2.0`, `pan_ukb_ld/pack` and `pan_ukb_ld/work` would free about 185 GB.

`$WORK_DIR/pan_ukb_ldscores` is a dangling symlink to the old `/dev/shm` dir and can be removed. The A2.2 rerun wrapper `$WORK_DIR/run_panel_a2.sh` (parallel ranged download from Dropbox) is not in the repo. If single-connection Dropbox downloads stay this slow, it's worth folding into `10_build_1000g_panel.sh`.

### 9.4 Outputs to look at

- **A1.1 allele-direction JSONs:** `/local-scratch/data/pan_ukb_ld/phase0/compare_{EUR,EAS,CSA,AFR,AMR}_chr22.json`
- **A2.2 panel build log:** `/local-scratch/data/genehackman_grch38/panel.log`. The first attempt's log is in the overnight tmux scrollback only.
- **A2.3 LD scores:** the log wasn't kept on disk; its summary is in the A2.3 status above.
- **A2.4:** `/local-scratch/data/genehackman_grch38/validation_report.md` (plus `validation_report_{EUR,EAS}.md` and `validation_{EUR,EAS}.log`)
- **QTL (B2):** `qtl_1.2.0/<dataset>/GRCh38/summary.tsv` and `qtl_1.2.0/<dataset>/GRCh37/liftover_summary.tsv`; logs `qtl_from_gpm.log` and `lift_qtl.log`

## 10. Uncommitted changes (download script, 5 Nov)

The A2 changes were committed by Andrew as `0e4b2e1`. New since then:

| File | Change |
|---|---|
| `download_reference_data.sh` (new, repo root) | the user download script (§11) |
| `reference_data/grch38/50_assemble_bundle.sh` | no GRCh37 bundle; prints upload commands that keep md5s (no composite uploads) and leave out the GRCh37 QTL copies |
| `snakemake/util/common.smk` | `stop_with_message()`: plain, boxed message with no traceback. `check_reference_data()`: lists everything missing (1000G panel or Pan-UKB pack per ancestry, RSID marker map, dbSNP for `populate_rsid: true`, LD scores, QTL dataset, GRCh37 files) with one download command. `validate_ancestries()` explains MID is Pan-UKB only. `parse_pipeline_input(output_builds=...)` |
| `snakemake/{standardise_gwas,disease_progression,finemap,coloc,compare_gwases,qtl_mr}.smk` | each calls `check_reference_data()` (compare_gwases with its LD scores, qtl_mr with its dataset); qtl_mr needs `output.build: GRCh38` and now validates ancestries against its own `ld_reference` |
| `snakemake/util/constants.smk`, `R/constants.R`, `R/data_conversions.R`, `R/liftover.R` | error messages and comments name `download_reference_data.sh` (`--grch37`, `--qtl`, dbSNP); `qtl_dataset_dir()` says the QTL data are GRCh38 only |
| `tests/testthat/test_reference_builds.R` | expects the `--grch37` message; skipped when the GRCh37 panel is installed (it failed against the live dir before) |
| `README.md` | §3 rewritten around the script: options, sizes, notes, Pan-UKB attribution; `.env` table |
| `CONTRIBUTING.md`, `tests/e2e_tests/run_test_pipelines.sh` | e2e test data: `--ancestries EUR,EAS --grch37 --qtl all` (~46 GB) |
| `PIPELINES.md`, `snakemake/{coloc,compare_gwases,disease_progression,finemap,qtl_mr,standardise_gwas}.md` | `output.build: GRCh37` needs `--grch37`; `qtl_mr` is GRCh38 only; removed the unimplemented `PIPELINE_DATA_DIR/qtl_datasets/` fallback; `finemap.md` panel path |
| `reference_data/grch38/{RUNBOOK,README}.md` | §4 GRCh37 QTL copy not published; §7 how the script lists the bucket, upload commands, Pan-UKB decision, testing with `file://`; §9 e2e data via the download script |
| `REVISION1_PLAN.md` | this update |

## 11. Download script (5 Nov 2026)

**`./download_reference_data.sh`** (bash + curl; the buckets are public over HTTPS, so no gsutil):

| Option | Default |
|---|---|
| `--dir DIR` | `PIPELINE_DATA_DIR` from `.env` |
| `--ancestries` | `EUR` (EUR EAS AFR AMR SAS MID, or `all`; SAS also fetches Pan-UKB CSA) |
| `--qtl`, `--qtl-dir` | `none`; `QTL_DATA_DIR` from `.env`, else `DIR/qtl_data` |
| `--no-dbsnp` | dbSNP included |
| `--grch37` | off; fetches the GRCh37 panel, LD scores (and dbSNP b37) from `1.1.0/` |
| `--dry-run`, `--verify`, `--jobs` (8), `--version` (1.2.0) | |

**How it works (no manifest).**
- It lists the bucket through the public GCS JSON API (`storage/v1/b/<bucket>/o?prefix=…`, paged). That gives each object's name, size and md5; the md5 is base64, decoded in awk.
- It picks files by path:

  | Path | Downloaded |
  |---|---|
  | `1000genomes/b38_dbsnp156/<POP>.{bed,bim,fam,frq}` | per ancestry |
  | `pan_ukb_ld/<POP>/` | per ancestry |
  | `LDSCORE/pan_ukb/<POP>/` | per ancestry |
  | `dbsnp/b38_dbsnp156/` | unless `--no-dbsnp` |
  | `liftover/`, `gene_name_map.tsv`, `marker_to_rsid.tsv.gz`, `pan_ukb_ld/MANIFEST.json`, `w_hm3.snplist` | always |
  | `1.1.0/…b37…` | with `--grch37` |
  | `<dataset>/GRCh38/` (QTL bucket) | with `--qtl` |
- It skips files already present with the right size, and checks free space where each part goes.
- It downloads in parallel to `.part` files with resume, checks size and md5, then moves each file into place.
- It checks that every selected file is in place, and prints the `.env` lines and what the data supports.
- **Testing before the upload:** `GENEHACKMAN_DATA_URL=file:///dir` lists with `find`, and checks sizes only.
- **Uploads must keep md5s.** Composite uploads have none, so `50_assemble_bundle.sh` prints `gsutil` commands with `parallel_composite_upload_threshold=0`.

**Sizes.** Every install includes 0.2 GB of core files.

| Item | Size (GB) |
|---|---|
| EUR | 7.8 |
| EAS | 7.3 |
| SAS | 8.7 |
| AMR | 16.9 |
| MID | 17.5 |
| AFR | 27.1 |
| dbSNP b38 | 11.7 |
| eqtlgen | 1.8 |
| brainmeta | 2.3 |
| ukb_pqtl | 0.4 |
| `--grch37` | ~1.9 per ancestry, + 11.4 for dbSNP b37 |
| everything | 124 |

A default EUR install is 19.7 GB.

**Tested (5 Nov), after removing the manifest.**
- **Bucket listing over HTTPS:** all 944 objects of `1.1.0/` listed. Size and md5 agree with an independent Python decode for 944/944, and the decoded md5s match local files. Paging gives the same result with 100-object pages (10 pages).
- **Classification:** the same files and sizes per group and ancestry as the earlier manifest (GRCh37: 21 panel, 25 dbSNP, 441 LD score files; full 1.2.0: 101.6 GB plus GRCh37).
- **Real download from a `file://` stand-in for 1.2.0** (symlinks to the live data): EUR plus `ukb_pqtl`, 20.1 GB, 2,475 files in 17 s. Byte-identical to the sources, and a rerun fetches nothing.
- **Real HTTPS download with md5 checks** (`--version 1.1.0`, core files). A same-size corrupted file is found and repaired by `--verify`. An unknown release gives a clear error.
- **Bugs found and fixed while testing:**
  - the awk join selected nothing when the target dir was empty;
  - a `--verify` message was written into the file list;
  - `find` didn't follow symlinks;
  - an empty md5 field shifted the `xargs` arguments (now `-` = size check only).
- **Pipeline messages** (Snakemake dry runs with the host `genehackman` conda env, Snakemake 7.26; the develop container's Snakemake fails on `distutils` under Python 3.12). Each case prints one boxed message, with no traceback, listing what is missing and one command:
  - EAS panel missing;
  - EAS pack missing;
  - GRCh37 panel and marker map missing (`--grch37`, plus `--no-dbsnp` because that dir has no dbSNP);
  - compare_gwases with EAS panel and LD scores missing;
  - qtl_mr with EAS panel and eqtlgen both missing (one command with `--ancestries EAS --qtl eqtlgen --qtl-dir …`);
  - `QTL_DATA_DIR` unset;
  - dbSNP missing for `populate_rsid: true`;
  - an empty data dir.

  qtl_mr with GRCh37 output, and MID with `ld_reference: 1000genomes`, also get plain explanations. All 19 e2e inputs still plan against the live `pipeline_data`.
- **R tests:** `test_reference_builds`, `test_data_conversions`, `test_mr` and `test_liftover` pass.
- **Not tested:** the published 1.2.0 data (nothing is uploaded yet), and macOS (the script falls back to BSD `stat`/`md5`).

**Test files** in `/local-scratch/data/genehackman_grch38/download_test/` (~20 GB: `data/`, `qtl/`, `nodbsnp/`, `empty/`, `https/`) can be deleted. The `file://` stand-in bucket is symlinks in the session scratchpad.

## 12. Pan-UKB redistribution: decision (5 Nov 2026)

**Decision (Andrew Elmore, 5 Nov 2026):** GeneHackman re-hosts the Pan-UKB-derived LD pack (`genomic_data/pan_ukb_ld/`) and the Pan-UKB LD scores (`LDSCORE/pan_ukb/`) in `gs://genehackman/1.2.0/`. The basis is the Pan-UKB downloads page's statement that its data are released openly.

**What the Pan-UKB pages say** (read on 5 Nov 2026):
- [Pan-UKB downloads page](https://pan.ukbb.broadinstitute.org/downloads):
  - Openness: *"All data here are released openly and publicly for the benefit of the wider biomedical community."* **This is the statement the decision rests on.**
  - The licence: *"The GWAS results data produced by the Pan-UKB are available free of restrictions under the Creative Commons Attribution 4.0 International (CC BY 4.0)."* This explicitly names the GWAS results only. LD matrices and LD scores are covered by the "released openly" statement above.
  - Terms: *"use of this data is bound by all terms of usage of the UK Biobank"*, and *"All users of this data agree to not attempt to reidentify participants."*
  - Attribution: *"We request that you acknowledge and give attribution to both the Pan-UKB project and UK Biobank, and link back to the relevant page, wherever possible."* and that publications cite Karczewski, Gupta, Kanai et al., *Nat Genet* 57, 2408–2417 (2025), [doi:10.1038/s41588-025-02335-7](https://doi.org/10.1038/s41588-025-02335-7).
  - Publication: *"From the perspective of the Pan-UKB team, there are absolutely no restrictions or embargoes on the publication of results derived from this data."*
- [Pan-UKB LD scores and matrices page](https://pan.ukbb.broadinstitute.org/docs/ld): no separate licence or terms.

**What GeneHackman does about the conditions:**
- The README's "Download reference data" notes say the pack and LD scores come from Pan-UKB, that UK Biobank terms apply, that users mustn't attempt re-identification, and that they should acknowledge Pan-UKB and UK Biobank and cite Karczewski et al. 2025.
- RUNBOOK §7 records the decision for maintainers.
- For the manuscript: credit Pan-UKB and UK Biobank and cite Karczewski et al. 2025 in the data-availability statement and the methods, and say that the LD pack is derived from the Pan-UKB LD matrices (restricted to ±1 Mb, eigenvalues clipped).

**Not covered by this decision:** the QTL bucket (UKB-PPP, eQTLGen, BrainMeta summary statistics, via the genotype-phenotype map). Their terms haven't been checked.
