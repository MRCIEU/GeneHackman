# GeneHackman revision 1: handoff

**State on 4 November 2026.** Written so that work can carry on on **ieup1**, where the data now lives: the Pan-UKB LD pack is 70 GB and doesn't fit on the laptop. It records what has been done, what is left in order, how to test, and the decisions still open.

The detailed planning documents are in the claude.ai project "GeneHackman" (`claude/GeneHackman_Revision1_Work_Plan.md`, `claude/revision1_status.md`, `claude/pan_ukb_ld_plan.md`, `claude/grch38_first_plan.md`, `claude/grch38_runbook.md`). Everything needed to continue is repeated here.

---

## 1. Where things are

### Code

All work is on **`feature/1.1.1-reviewer-feedback`**, on top of `main` (`2682858`). Push from the laptop, pull on ieup1.

| Commit | What |
|---|---|
| `b9c6247` | Audit correctness fixes: B1 opengwas preset, R1-M3 build-aware rsIDs, B22, B23, B10/B11, B16, B17, B8 (see §5) |
| `bbc675a` | Pan-UKB LD pack: build scripts (`reference_data/pan_ukb_ld/`), readers (`R/ld_reference.R`, `scripts/ld_reference.py`), LD-matrix clumper (`R/clumping.R`), `ld_reference: pan_ukb` in every LD-using pipeline |
| `51fdcba` | Pack: variants sharing a GRCh37 key (e.g. an insertion and a deletion `A/AT` at one position) get no keys |
| `bbcf276` | Pack radius ±1 Mb; settings-aware chromosome markers; PSD fix by eigenvalue clipping instead of a ridge |
| `c838069` | End-to-end tests run every pipeline with `ld_reference: pan_ukb` too |
| `0e6b977` | GRCh38-first: liftover QC, build-aware reference data, default output build GRCh38, `reference_data/grch38/` scripts and `RUNBOOK.md` |
| `aaede5c` | Pack build: Hail no longer retries forever when its Java process dies; extraction restarts itself from checkpoints |
| (this file) | Handoff notes |

### Data on ieup1

| Path | What |
|---|---|
| `/local-scratch/data/pan_ukb_ld/` | Pack build `WORK_DIR`: `pack/` (the result), `raw/`, `state/` (stage markers), `logs/`, `env/` (Hail environment), `phase0/` |
| `/local-scratch/data/pan_ukb_ld/pack/` | **Finished pack, 70 GB.** Sizes: EUR 5.6 GB, EAS 5.3, CSA 6.4, AMR 15, MID 17, AFR 22 |
| `/local-scratch/data/pan_ukb_ld/genehackman_develop.sif` | GeneHackman container (used as `R_CMD` for the R build steps) |
| `/home/wt23152/GeneHackman` | GeneHackman checkout the pack was built from |
| `/local-scratch/data/genehackman_grch38/` | Planned `WORK_DIR` for the GRCh38 reference data (`reference_data/grch38/config.sh`); not created yet |

The Pan-UKB source blocks (downloaded LD) were deleted after each ancestry (`KEEP_SOURCE_BLOCKS=0`). Anything that changes the LD itself (radius, dead-zone, quantisation) would mean downloading again. The GRCh38 keys can be redone without downloading (§3, step A3).

---

## 2. What has been done

### 2.1 Audit correctness fixes (`b9c6247`)

| ID | Fix |
|---|---|
| B1 | `opengwas` preset: EA = ALT, OA = REF (ES and AF refer to ALT); signs were flipped |
| R1-M3 | RSIDs are looked up in the output build, and the `chr:pos` parser works row by row |
| B22 | Expected-vs-observed uses a two-sided observed p |
| B23 | Comparisons no longer overwrite each other with 4 or more GWAS |
| B10/B11 | SuSiE writes PIPs; overlapping windows are merged into regions; region and credible-set tables |
| B16 | coloc no longer zero-fills missing log Bayes factors |
| B17 | QTL coloc uses the QTL's own `study_type`, N and priors |
| B8 | Default clumping `--clump-r2 0.01 --clump-kb 1000` when not given |

### 2.2 Pan-UKB LD reference (`ld_reference: pan_ukb`)

**Source.** Pan-UKB LD matrices: Hail BlockMatrix, GRCh37, ±10 Mb, signed r on the ALT allele, 43 TB, from `s3://pan-ukb-us-east-1/ld_release/`. Sample sizes per ancestry:

| EUR | CSA | AFR | EAS | MID | AMR |
|---|---|---|---|---|---|
| 420,542 | 8,876 | 6,637 | 2,709 | 1,599 | 982 |

**Pack format ("GeneHackman LD pack v1").** Per ancestry, there is `info.tsv` plus, per chromosome, `variants_chr<C>.tsv.gz` and `ld_chr<C>.bin` with its `.chunks.tsv`.

- **Variants kept:** MAF ≥ 1% in that ancestry.
- **Pairs stored:** within ±1 Mb, as a ragged upper-triangle band.
- **Values:** r × 127 as int8; |r| < 0.02 is stored as 0; xz-compressed in chunks of 2,048 rows.
- **Orientation:** r is oriented to the alphabetically smaller allele (GeneHackman's EA), so it can be used directly with standardised GWAS z-scores.
- **Keys:** `SNP37` and `SNP38`, in GeneHackman's `CHR:BP_EA_OA` form. GRCh38 keys come from GeneHackman's own liftover, so a lifted GWAS matches exactly.
- **`LIFT_STATUS`** says why a variant has no GRCh38 key. The full description is in `reference_data/pan_ukb_ld/PACK_README.md`.

**Build.**

- `reference_data/pan_ukb_ld/run_in_tmux.sh phase0`, then `full`.
- Downloads only the LD blocks it needs, one ancestry at a time.
- Resumable: checkpoints per chunk, and stage markers that record the settings.
- Each chunk is validated byte for byte against the source (`state/<POP>.validate.json`).

**Phase 0 findings** (at ±500 kb):

- Every ancestry's chr22 was byte-identical to the source.
- Download was about 5.5 TB at ±500 kb (more at ±1 Mb).
- In the CDKN1A test region, 1.6% of the LD r² mass lies beyond 500 kb, including pairs with |r| up to 0.8. This is why the radius became ±1 Mb.

**Size.**

- The final pack is 70 GB.
- Lossless tricks (stronger xz, bzip2, delta coding) gain nothing; tested on 1000G-derived LD run through the same encoder.
- Lossy options (coarser rounding, a larger dead-zone) save 30–57% but changed SuSiE results in simulations, so they were not used.
- The noise-dominated ancestries (AMR, MID, AFR) are most of the size.

**Pipelines.** With `ld_reference: pan_ukb`, every LD step uses the pack:

- **Clumping:** `R/clumping.R`, which follows PLINK 1.9 `--clump` rules on the stored r. PLINK parity is tested.
- **SuSiE, MultiSuSiE, LocusZoom r² and EAF filling.**
- **Matching** is on SNP keys, not RSIDs.
- **Regions table:** `finemap_regions.tsv` records `LD_REFERENCE`, `N_REF`, `FRAC_BEYOND_LD_RADIUS` and `LD_MIN_EIGENVALUE`.

**Bug found and fixed (`bbcf276`): positive semi-definite (PSD) correction.**

- Rounded LD has negative eigenvalues: −0.1 to −1.0 for 1 Mb regions in 1000G-based tests.
- The old fix added |λmin| to the diagonal. In simulations this lost or invented credible sets (sets matched exact LD in only 3–53% of runs).
- The new fix clips negative eigenvalues and rescales to a unit diagonal. Results are close to exact LD, with causal-variant coverage unchanged.

**Clumping caveat for the methods.** PLINK 1.9 clumps on haplotype (EM) r², while Pan-UKB stores genotype-dosage r. Results therefore differ slightly for pairs at the threshold, even before the panels differ.

### 2.3 Liftover QC (work plan §3.9, `0e6b977`)

**Measured on the old liftover, 1000G EUR (8.55M variants).**

- Silently wrong:
  - 19,418 minus-strand variants were placed 2 bp too low, with alleles not complemented;
  - 2,498 conversion-unstable positions were kept.
- Dropped without a reason: 269 chromosome changes and 2,765 alt-contig mappings.

**Now, in `R/liftover.R`.**

- **Positions:** correct one-base BED intervals, joined back on a row ID.
- **Dropped variants are counted by reason:** `unmapped`, `non_primary_contig`, `chrom_change`, `strand_flip` (or complemented if non-palindromic, as an option), `round_trip` (conversion-unstable) and `duplicate_target`.
- **Report:** written beside each standardised GWAS. The run stops if more than 5% is dropped.
- **Declared-build check:** before lifting, positions are checked against the pack's GRCh37 and GRCh38 keys.
- **GRCh36 to GRCh38** is supported directly (hg18ToHg38).
- **Tests:** unit tests use the real UCSC chains (set `GH_LIFTOVER_CHAIN_DIR`).

The hg38 test fixtures were regenerated, because the old ones contained the mis-lifted variants.

### 2.4 GRCh38-first (`0e6b977`)

- **Default output build is GRCh38.** The input build still defaults to GRCh37. `output.build: GRCh37` is opt-in and needs the GRCh37 extras.
- **Every reference resource is chosen by output build:**
  - the 1000G panel, `.frq` and marker map (`genomic_data/1000genomes/b38_dbsnp156` or `b37_dbsnp156`);
  - dbSNP;
  - the QTL data (`<dataset>/<build>/...`);
  - the LocusZoom gene tracks (EnsDb v86 or v75).
- **Clumped positions** are taken from the GWAS by RSID, so they never depend on the panel's build.
- **LDSC:** `ldsc.ld_scores` is `pan_ukb` (default) or `1000genomes_b37`. LDSC matches on rsID, so either works with either build.
- **Pre-standardised `*_std.tsv.gz` inputs** must declare `build:` equal to `output.build`.
- **Pack keys:** the pack annotation now uses the package liftover. `run_build.sh update-keys` refreshes a finished pack's GRCh38 keys without touching its LD.
- **Data scripts** are in `reference_data/grch38/`: the 1000G 30x GRCh38 panel, Pan-UKB LD scores, QTL liftover, validation report and bundle assembly. **Instructions: `reference_data/grch38/RUNBOOK.md`.**

### 2.5 Verification so far (cloud, 1–3 Nov)

- testthat: 361 expectations pass. Python unittest: 9 pass.
- Synthetic Pan-UKB-format build passes end to end, including:
  - an interrupted and resumed run;
  - a phase 0 at another radius being rebuilt;
  - `update-keys`.
- Key parity: pack keys equal `standardise_gwas()` keys in both builds (74,953 of 74,953 both ways).
- PLINK clump parity on the same LD.
- Snakemake dry runs of all 17 end-to-end inputs, plus the error cases.
- R CMD check: no new notes.
- **Not yet run anywhere:** real end-to-end pipelines with the real pack or the GRCh38 data.

---

## 3. What is left, in order

### A. Finish and check the Pan-UKB pack (ieup1)

1. **Check the build finished and validated:**
   ```bash
   cd /home/wt23152/GeneHackman
   bash reference_data/pan_ukb_ld/run_build.sh status
   for p in EUR AFR AMR CSA EAS MID; do python3 -c "import json;d=json.load(open('/local-scratch/data/pan_ukb_ld/state/$p.validate.json'));print('$p', d['ok'], sum(len(c['problems']) for c in d['chromosomes'].values()))"; done
   bash reference_data/pan_ukb_ld/run_build.sh finalise    # MANIFEST.json with checksums
   ```
2. **Allele-direction check against 1000G. Not run yet, and important.**
   - Pan-UKB doesn't document which allele r and AF refer to. Validation only proved the copy matches the source.
   - Expect EAF r ≈ 1, `share_closer_if_swapped` ≈ 0, and sign concordance ≥ 0.98 for |r| > 0.3.
   - Needs the GRCh37 1000G files (`genomic_data/1000genomes/b37_dbsnp156/{EUR,EAS,SAS,AFR,AMR}.*`). `R_CMD`'s container has `plink1.9`; wrap it as in the earlier instructions.
   ```bash
   source reference_data/pan_ukb_ld/config.sh
   G1K=/path/to/1000genomes/b37_dbsnp156
   SIF=/local-scratch/data/pan_ukb_ld/genehackman_develop.sif
   PLINK_WRAPPER=$WORK_DIR/plink1.9
   printf '#!/usr/bin/env bash\nexec apptainer exec --bind /local-scratch,%s %s plink1.9 "$@"\n' "$G1K" "$SIF" > $PLINK_WRAPPER
   chmod +x $PLINK_WRAPPER
   for pair in EUR:EUR EAS:EAS CSA:SAS AFR:AFR AMR:AMR; do pop=${pair%%:*}; g=${pair##*:}
     $WORK_DIR/env/bin/python reference_data/pan_ukb_ld/61_compare_1000g.py --pack $WORK_DIR/pack --pop $pop --chrom 22 \
       --bfile $G1K/$g --frq $G1K/$g.frq --window 20000000-25000000 --plink "$PLINK_WRAPPER" \
       --out $WORK_DIR/phase0/compare_${pop}_chr22.json; done
   ```
3. **Refresh the GRCh38 keys with the new liftover.** This adds the `round_trip` and `non_primary_contig` statuses, and runs on code `0e6b977` or later. LD is untouched.
   ```bash
   bash reference_data/pan_ukb_ld/01_setup_env.sh          # also fetches hg38ToHg19 for the round trip
   bash reference_data/pan_ukb_ld/run_build.sh update-keys # all ancestries
   ```
4. **Optional:** re-run the truncation fine-mapping check (`62_truncation_finemap.R`) at ±1 Mb on a HERMES locus.
5. **Optional, pack size.** Measure lossless or near-lossless re-encodings on the real pack, decoding `ld_chr*.bin` and re-encoding, with no download. No script exists yet.
   - Section 2.2 says the gains are expected to be small.
   - **Per-ancestry download** is the practical answer to the size: `gsutil rsync` only `pan_ukb_ld/EUR` (5.6 GB). Document it in the README.

### B. GRCh38 reference data: follow `reference_data/grch38/RUNBOOK.md`

Steps 0–9 in the runbook:

1. second checkout and settings;
2. check the dbSNP b38 `.fst` files (Andrew has them; locate them first);
3. 1000G 30x GRCh38 panel;
4. Pan-UKB LD scores;
5. QTL data in GRCh38;
6. pack `update-keys` (same as A3);
7. validation report;
8. assemble bundle 1.2.0 and the GRCh37 extras;
9. LDSC benchmark;
10. end-to-end tests.

The runbook gives the expected output at each step.

**Unverified from Claude's side:** the download URLs, the Pan-UKB LD-score archive layout and MetaBrain's build. The scripts check and print these.

### C. Decisions needed (Andrew)

1. **LD reference coherence:** see §4 below. The recommendation is to make each run self-contained in one reference, including LD scores and RSID annotation.
2. **Default `ld_reference`.** It is still `1000genomes`. Flipping it to `pan_ukb` waits for the licence check: Pan-UKB is CC BY 4.0, "usage may be restricted by UK Biobank". Redistributing the pack from `gs://genehackman` needs confirming.
3. **Whether to keep building the 1000G GRCh38 panel** if Pan-UKB becomes the default. See §4: it is still the route for in-sample LD and a small download.

### D. End-to-end tests and release

`tests/e2e_tests/run_test_pipelines.sh`:

- Runs every pipeline with both LD references, in GRCh38, plus the GRCh37 opt-in runs.
- Needs bundle 1.2.0, the GRCh37 extras and the GRCh38 QTL data (runbook step 9).
- The container needs `EnsDb.Hsapiens.v86`, so rebuild it: `docker/requirements.R` is updated.

For the release:

1. Bump `DESCRIPTION` to 1.2.0. The image tag and bucket prefix follow it.
2. Push the image.
3. Upload the bundles.
4. Commit `tests/testing_complete.txt`.

### E. Remaining code work from the revision 1 work plan

Not started unless marked. IDs are from the work plan.

**LD and fine-mapping (R1-M2, m11, m12)**

- **C2.1:** user-supplied in-sample LD (bfile/pgen or precomputed matrices), the preferred tier. The provider dispatch in `get_ld_matrix()` is where it plugs in.
- **C2.3:** `susieR::kriging_rss` LD-mismatch QC, with a flipped-allele positive control.
- **C2.4:** L = 1 when N / n_ref > 10, using N_eff for binary traits. `N_REF` is already recorded.
- **C2.5:** use the median per-SNP N, not `N[1]` (`R/finemap.R`).
- **C9.1:** fine-mapping HTML report and supplementary credible-set table. PIPs and credible-set summaries are done (B10).
- **C9.2 / m12:** one region module for SuSiE and MultiSuSiE (`merge | none | tile`, `max_region_mb`). MultiSuSiE has no region cap yet.
- **B13:** check the CDKN1A five-credible-set pattern for LD mismatch: L = 1, kriging, 1000G vs Pan-UKB.
- **B14:** flag MHC and long-range-LD regions.
- **B15:** one documented `window_kb` default.

**Harmonisation (R2-1)**

- **B2:** strand and palindrome handling against reference alleles. The pack's EAF can serve as the reference AF.
- **B3:** the Z-only input path.
- **B6:** per-GWAS QC report.
- **B7:** P = 0, and D/I indels.
- **C3.3:** GRCh38 regression test plus MungeSumstats concordance.
- **C3.4:** a `gwas_catalog` (GWAS-SSF) preset.
- **§5:** known-answer harmonisation tests in CI.

**LDSC (R1-M1, m1, R2-2)**

- **C1.2:** benchmark against the Alkes Price group's scores (runbook step 8).
- **C5.1:** all pairs. `--rg a,b,c` only does the first against the rest.
- **C5.2:** liability-scale output with N_eff (`--samp-prev` / `--pop-prev`).
- **C5.3:** parsed LDSC table with flags for h² z < 4 and |rg| > 1.

**Binary traits (m4, m5, m7)**

- **C7.1–C7.4:** per-GWAS `study_type`, `N_case`, `N_control`, `prevalence` and N_eff; per-SNP N takes precedence; legacy coloc uses `type = "cc"`.

**Other audit items**

| Item | Fix |
|---|---|
| B9 | `|| echo headers` in `clumping_rule.smk` hides PLINK failures |
| B27 | SlopeHunter default never applies (AttributeError) and the manuscript says SlopeHunter; the demo used CWLS |
| B28–B31 | Collider-bias method details, labels and `throw()` |
| B18 | coloc prior sensitivity |
| B20 | Steiger test, and report MR hits that were not coloc-tested |
| B24–B26 | Winner's curse, Q-test caveats, tidy-ups |
| B32 | Pin package versions and cite the image digest |
| B34 | e2e tests in CI with output assertions |
| B35 | Telemetry: make it opt-in and document it |
| B36–B38 | Count skipped loci; DESCRIPTION placeholder; munge `sed` patch |

Also not yet done: both-direction comparisons, and the `subset(mr_results, exposure %in% ...)` bug when `exposures` is set.

**Liftover L8 (later):** a reference-allele check against the target FASTA, or BCFtools/liftover as an engine.

### F. Re-analysis, manuscript and response letter

- **Work plan §6:** the HF demonstration re-run, including:
  - correct N and N_eff, and `study_type`;
  - strict clumping;
  - LDSC all pairs with the new LD scores;
  - fine-mapping with L = 1 vs 10, kriging, and 1000G vs Pan-UKB;
  - coloc after the B16/B11 fixes;
  - disease progression with SlopeHunter plus CWLS;
  - **E1** interpretation paragraphs for the multi-ancestry and disease-progression components.
- **Response letter:** the M1, M2, M4 and R2-1 sections predate the co-author comments and the Pan-UKB and GRCh38 decisions. The M3 text overstates the `chr:pos:a1:a2` parser: EA/OA must come from their own columns.
- **Methods sentences** to add:
  - Pan-UKB genotype r vs PLINK haplotype r²;
  - the PSD clipping;
  - liftover QC counts;
  - GRCh38 default.

---

## 4. LD references and ancestries: how they fit together

### One reference per run: what the code does now

`ld_reference` is a single setting for the whole run. Every step that uses LD reads from that one source, so a run never clumps with one reference and fine-maps with another:

| Step | `ld_reference: 1000genomes` | `ld_reference: pan_ukb` |
|---|---|---|
| Clumping | PLINK `--clump` on the 1000G panel | `R/clumping.R` on the pack (same rules) |
| SuSiE, MultiSuSiE | PLINK `--r` on the panel | pack |
| LocusZoom r² | panel | pack |
| EAF filling | panel `.frq` | pack `EAF` |

**Where the other source is still used today.** These steps don't use LD between variants:

- **LDSC LD scores:** `ldsc.ld_scores` defaults to `pan_ukb` whatever `ld_reference` is.
- **Partial RSID annotation:** uses the 1000G marker map whatever `ld_reference` is (RSIDs only).
- **Declared-build check:** uses the pack's keys (positions only).

**Proposal** (not implemented; needs Andrew's decision): make every run fully self-contained in its reference.

- `ldsc.ld_scores` would default to the scores of the chosen `ld_reference`. For `1000genomes`, the GRCh37-built scores are rsID-keyed and work for GRCh38; move them from the extras into the main bundle.
- With `pan_ukb`, partial RSIDs would come from the pack's `RSID_UKB` column. A `pan_ukb` run would then need no 1000G files at all.
- MultiSuSiE would keep using one reference for all ancestries in a run, which it already does.

### Ancestry codes

| Code | 1000 Genomes (2,504 unrelated; phase 3 for GRCh37, 30x for GRCh38) | Pan-UKB (UK Biobank participants grouped by genetic similarity) |
|---|---|---|
| EUR | 503: CEU, TSI, FIN, GBR, IBS | 420,542 |
| EAS | 504: CHB, JPT, CHS, CDX, KHV | 2,709 |
| AFR | 661: YRI, LWK, GWD, MSL, ESN, ASW, ACB | 6,637 |
| AMR | 347: MXL, PUR, CLM, PEL | 982 |
| SAS (1000G) / CSA (Pan-UKB) | 489: GIH, PJL, BEB, STU, ITU | 8,876 (GeneHackman maps SAS → CSA) |
| MID | — | 1,599 (Pan-UKB only) |

**What the labels mean.**

- The labels are continental, and the two resources sample them differently. 1000G samples named populations in their home countries. Pan-UKB groups UK residents by genetic similarity.
- The match between a GWAS and its LD panel is what matters, not the label. Examples:
  - for HERMES (mostly European cohorts), Pan-UKB EUR is far larger and well matched;
  - for BBJ (Japanese), 1000G EAS includes JPT, while Pan-UKB EAS is UK residents of mostly Chinese ancestry;
  - neither is in-sample, and the C2.4 L = 1 rule matters for both.
- AMR in both is admixed and small. Treat AMR fine-mapping as low-confidence with either reference.
- The pipeline validates codes against the chosen reference: MID only with `pan_ukb`, and SAS/CSA accepted for `pan_ukb`.

### Why keep both

- **1000G (PLINK path).**
  - It is the same code path that user-supplied in-sample bfiles (C2.1) will use. In-sample LD is the reviewers' preferred option.
  - It is a small download (about 10 GB).
- **Pan-UKB.**
  - The large panels the reviewers asked for. EUR is 836 times larger than 1000G EUR.
  - Six ancestries.
  - 70 GB in total, downloadable per ancestry.

---

## 5. How to test

```bash
# R unit tests (from the repo root; R with GeneHackman's dependencies, e.g. the container)
GH_LIFTOVER_CHAIN_DIR=/path/to/genomic_data/liftover \
  Rscript -e 'testthat::test_local(".", reporter = "summary")'     # liftOver must be on PATH for the real-chain tests

# Python tests
cd tests/python && python -m unittest -q && cd ../..

# Pan-UKB pack build on a synthetic source (needs the Hail env)
PYTHON=/local-scratch/data/pan_ukb_ld/env/bin/python R_CMD="$R_CMD" LIFTOVER_BIN=$(which liftOver) \
  CHAIN_FILE=/path/to/genomic_data/liftover/hg19ToHg38.over.chain.gz \
  bash reference_data/pan_ukb_ld/tests/run_synthetic_build.sh /tmp/synthetic_pack_test

# Key parity: pack keys == standardise_gwas() keys in both builds
WORK_DIR=/tmp/key_parity GENEHACKMAN_DIR=$PWD CHAIN_FILE=.../hg19ToHg38.over.chain.gz LIFTOVER_BIN=$(which liftOver) \
  Rscript reference_data/pan_ukb_ld/tests/key_parity.R

# Snakemake dry run of one pipeline
snakemake -n -p --snakefile snakemake/finemap.smk --config genehackman_input=tests/testthat/data/snakemake_inputs/finemap_pan_ukb.yaml

# End-to-end (real data, container): see §3 D
./tests/e2e_tests/run_test_pipelines.sh
```

---

## 6. Gotchas

- **Don't `git pull` in a checkout while a bash build script from it is running.** Bash reads scripts as it goes. Use `git worktree add` for a second checkout.
- **Stage markers record their settings.** Changing `RADIUS_KB`, `DEADZONE`, `MAF_MIN` or `CHUNK_ROWS` rebuilds LD chromosomes (from new downloads). `FETCH_ROWS`, `HAIL_DRIVER_MEMORY` and `EXTRACT_ATTEMPTS` can change freely.
- **Hail's Java process stopped once** with a normal shutdown, on 1 Nov at 14:29. Hail then retried "Connection refused" for two days.
  - Since `aaede5c`, Hail gives up after 5 tries, and the build restarts the extraction itself from the checkpoints.
  - The cause of that shutdown wasn't confirmed. It wasn't out of memory: the machine has 1 TB of RAM.
- **R steps in the build** run through `R_CMD` (the apptainer container) with `pkgload::load_all()` of the checkout. The container only provides R packages, so new R code needs no rebuild. The exception is new packages such as EnsDb v86.
- **Liftover needs both directions' chains** (`hg19ToHg38` and `hg38ToHg19`) for the round-trip check, plus `hg18ToHg38` for GRCh36 input.
- **A "declared as GRCh37, but its positions match GRCh38" stop when standardising** means the GWAS's `build:` is wrong. That's the check working.
