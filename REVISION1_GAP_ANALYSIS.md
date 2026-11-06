# Revision 1: reviewer points, status and issues the reviewers didn't raise

6 November 2026. Code state: branch `feature/1.1.1-reviewer-feedback` (local, not pushed), commits on top of `168f299`:

| Commit | Content |
|---|---|
| `ac0f629` | Binary traits and N_eff; per-variant N takes precedence; Z-only input. Fine-mapping: user/in-sample LD, allele alignment, LD-mismatch QC, automatic L, region strategies, MHC flag (SuSiE and MultiSuSiE). |
| `03001e5` | LDSC: every pair within ancestry, effective N, liability-scale h², parsed tables with reliability flags. |
| `a9403af` | LD-mismatch flags per credible set and in the HTML reports; fine-mapping HTML report; LD-mismatch QC in MultiSuSiE (Python port of susieR's diagnostics, tested against susieR); MungeSumstats script with Reviewer 1's chr:pos case. |
| `941dbc9` | Strand/palindrome check against the LD reference; known-answer harmonisation tests; MungeSumstats concordance script; case-control coloc datasets; per-GWAS labels; B9, B27, B35, B32 (versions in reports); documentation. |

Tests: testthat 138 blocks pass (110 before; +28), Python 17 pass with MultiSuSiE installed, Snakemake dry runs of finemap (single and multi-ancestry), compare_gwases, coloc and disease_progression parse with the new options.

**Status key.**

- **Done:** code, tests and docs are in.
- **Code done, needs re-analysis:** the letter's numbers come from the HF re-run on ieup1.
- **Partial:** part of what the draft letter promised is in; the rest needs a decision, or the letter changed.
- **Open:** nothing yet.

---

## 1. Every reviewer point

### Reviewer 1

| Point | Asked | Status | Evidence | Still to do |
|---|---|---|---|---|
| Overall: positioning | Place GeneHackman among harmonisation and post-GWAS tools | **Done (manuscript, short)** | One Background paragraph: dedicated tools (MungeSumstats, GWAS-VCF) standardise as an end in themselves; GeneHackman's standardisation is the first step of its own downstream analyses | Andrew (6 Nov): no tool-by-tool comparison or supplementary table. Gibran's broader suggestion (MR pipelines, vl.genoma.io, FUMA) not taken up; check he agrees. |
| Overall: best practice | Fine-mapping and LDSC to follow best practice | **Done** | See M1, M2, m1, m4, m7, m12 | – |
| Overall: critical evaluation | Results presented without critical evaluation | **Open (re-analysis + text)** | Placeholders in the manuscript | Gibran: "high level interpretation of each component" (work plan E1). |
| M1 | LD scores built with `--ld-wind-kb 1`; benchmark against public scores | **Code done, needs benchmark** | Pan-UKB LD scores are the default (`ldsc.ld_scores: pan_ukb`, 1.2.0 bundle); `reference_data/grch38/20_import_pan_ukb_ldscores.py` | Runbook step 8: benchmark against the Alkes Price EUR/EAS scores (needs HF GWAS paths). Gibran asked to show that we reproduce existing scores. The draft's "ldsc.ref_ld_chr / w_ld_chr" YAML options don't exist: the letter now says `ldsc.ld_scores`. Confirm the `--ld-wind-kb` flag in the v1.1.0 `.log` files. |
| M2 | In-sample LD option, or explain; state that 1000G LD is used | **Done** | `ld_bfile:` per GWAS (PLINK), used by SuSiE and MultiSuSiE, matched on position and alleles; `max_causal: auto` (L = 1 when N > 10 × N_REF); `ld_mismatch_qc()` (estimate_s_rss + kriging_rss); allele alignment for 1000G; Pan-UKB pack as a large precomputed public LD; `finemap.md` "LD reference and its limits" | Precomputed **user** matrices (.bcor/.npz) aren't supported: the letter says PLINK files plus the Pan-UKB pack. gnomAD LD (Andrew's and Gibran's comments) is replaced by the Pan-UKB pack. Re-run HF with L = 1 vs 10, kriging, 1000G vs Pan-UKB (B13: CDKN1A). |
| M3 | Positioning; chr:pos:a1:a2 to GRCh38 + rsID failed | **Done (code), Partial (letter)** | rsIDs looked up in the output build (b9c6247); the `chr:pos` parser fills CHR/BP row by row; known-answer test | The parser does **not** take alleles from the identifier: EA/OA must come from columns (the draft said it does; corrected). Gibran: "check if it's a bug?" Reproduce the reviewer's case if possible. MungeSumstats concordance: run `scripts/validation/compare_with_mungesumstats.R` on ieup1. |
| M4 | Document how to add modules (CARMA, SLALOM) | **Done (docs)** | `CONTRIBUTING.md` "Adding a method to an existing pipeline": stage interfaces + CARMA worked example | No `finemap.method` option yet: the letter describes it as the route, not a shipped feature. Gibran: make it more prominent (Discussion sentence). |
| m1 | LDSC red flags; all pairs; prevalence; compare with Henry et al.; consider LDAK | **Code done, needs re-analysis** | `run_ldsc.sh` (all pairs, N_EFF, --samp-prev/--pop-prev), `ldsc_results_tables()` (h² z < 4, \|rg\| > 1), report tables | Re-run; comparison table with Henry et al. 2025. LDAK isn't added (letter: listed as future module). |
| m2 | Adoption in biobank TREs | **Done (docs), letter to confirm** | `finemap.md` "In-sample and user-supplied LD" | Confirm container support on UKB-RAP / All of Us / FinnGen before naming them. |
| m3 | GRCh38 | **Done** | GRCh38 default output; 1000G 30x GRCh38 panel; dbSNP b38; EnsDb v86; QTL in GRCh38; liftover QC | Container rebuild (EnsDb v86); upload 1.2.0. |
| m4 | Binary traits throughout | **Done** | `study_type`, `N_case`, `N_control`, `prevalence`; N_EFF in standardisation; LDSC liability scale; SuSiE/MultiSuSiE use N_EFF; case-control `coloc_dataset()`; QTL coloc uses QTL's own type (B17) | coloc.susie (the main coloc path) uses fine-mapping output, so trait type enters through N_EFF only; the letter says so. |
| m5 | N = 1.2M placeholder in coloc example | **Open (re-analysis)** | – | Re-run with N_EFF; regenerate Fig 2d; config fix in Supplementary. |
| m6 | Acknowledge sample overlap | **Done (docs/report), Open (manuscript)** | Static warning in coloc report; `coloc.md`; LDSC `GCOV_INTERCEPT` | The draft's "report warns when LDSC intercept indicates overlap / user declares overlap" isn't implemented; the letter now describes the static warning. |
| m7 | Effective N; sum of cohort N_eff | **Done** | `populate_effective_n()`: N_EFF column kept, else N_CASES/N_CONTROLS, else N × 4φ(1−φ), else study counts | Confirm whether HERMES files carry per-variant N_eff. |
| m8 | "multi-trait colocalization" | **Done (manuscript)** | Tracked change | HyPrColoc/moloc: mentioned as extensions, not implemented. |
| m9 | BBJ not in GWAS comparison | **Open (decision)** | – | Option A (add) or B (correct text). Placeholder in manuscript. |
| m10 | Table 1 wording | **Done (manuscript)** | Tracked change | Parameter table (Supplementary). |
| m11 | Fine-mapping results not reported | **Done** | HTML report (`finemap.Rmd`, commit a9403af); PIPs, credible-set table with `LD_FLAG` (OUTLIER_IN_CS, LD_MISMATCH, OUTLIERS_IN_REGION, MHC); regions table with L, LD_MISMATCH_S, N_LD_OUTLIERS; coloc and QTL MR reports list flagged credible sets; coloc results have `ld_check` | Fill the numbers from the HF re-run. |
| m12 | Region definition: max size, no-merge, tiling | **Done** | `region_strategy: merge | none | tile`, `max_region_kb`, same rules in SuSiE and MultiSuSiE, MHC flag, documented | – |
| m13 | Duplicate reference | **Done (manuscript)** | Tracked change | – |
| m14 | Cite Chang et al. 2015 | **Done (manuscript)** | Tracked change | – |
| m15 | Shah 2020 vs Henry 2025 | **Done (manuscript)** | Tracked change | – |
| m16 | HERMES name | **Done (manuscript)** | Tracked change | – |
| m17 | Describe datasets; phenotype names not file names | **Code done, needs re-analysis** | `label:` per GWAS used in coloc, LDSC and heterogeneity outputs | Dataset table (placeholder); re-made figures. Miami/Manhattan plot titles still use file names. |
| m18 | Colocalization/colocalisation | **Done (manuscript)** | Tracked change; report title | – |
| m19 | Column auto-detection; BBJ columns | **Done** | Docs state no auto-detection; `gwas_catalog` preset; per-variant N precedence | Reconcile BBJ N (YAML 203,234 vs file 206,495; 16,251 cases / 197,577 controls). |
| m20 | Local file paths in configs | **Open (supplementary)** | Repo templates are clean | Clean the demo YAMLs in the Supplementary/Zenodo record. |

### Reviewer 2

| Point | Asked | Status | Evidence | Still to do |
|---|---|---|---|---|
| R2-1 | Show harmonisation is correct (effect direction, ambiguous variants) | **Done (tests), needs concordance run** | `test_harmonisation_known_answers.R` (12 cases: swapped alleles, flip_alleles, case, chr prefix/X, indels and I/D, OR, Z-only, LOG_P, P = 0, chr:pos IDs, per-variant N, presets, strand/palindrome); `harmonise_against_reference()` with `<output>.allele_qc.tsv`; liftover round-trip tests | Run the MungeSumstats script on HF and BBJ for the concordance numbers. TwoSampleMR comparison: not scripted (letter drops it). Known-answer tests run in the testthat suite (CI runs it). Gibran: use the R1 checks to answer this; done. |
| R2-2 | Which LDSC pairs and LD scores; EUR vs Japanese? | **Code done, needs re-analysis** | All pairs within ancestry; LD score set recorded | Popcorn not added: the letter states the limitation. |

---

## 2. Issues the reviewers did not raise

Found during the audit, the Pan-UKB/GRCh38 work and today's implementation. "Fixed" means code and tests are in.

| # | Issue | Effect on the submitted results | Status |
|---|---|---|---|
| 1 | **OpenGWAS VCF preset** took the effect allele as REF; effects and frequencies refer to ALT (B1) | Signs reversed for files read with this preset; not used in the manuscript | Fixed (b9c6247) |
| 2 | **Expected vs observed replication** counted observed hits with a one-sided p against a two-sided expectation (B22); with ≥ 4 GWAS some comparisons overwrote others (B23) | Observed replication inflated in the comparison pipeline | Fixed |
| 3 | **coloc** gave a log Bayes factor of 0 to variants missing from one dataset (B16) | Posterior probabilities biased where variant sets differ | Fixed |
| 4 | **QTL coloc** used the GWAS's N, trait type and priors instead of the QTL's (B17) | Wrong priors for QTL coloc | Fixed |
| 5 | **Clumping defaults:** with only `--clump-p1`, PLINK's r² 0.5 / 250 kb applied (B8); PLINK failures were turned into an empty result by `|| echo` (B9) | Several correlated leads per locus; failures looked like "no hits" | Fixed (B9 today) |
| 6 | **Overlapping single-ancestry fine-mapping windows** produced duplicate signals (B10/B11), and MultiSuSiE merged windows with no size limit | Duplicate entries in Fig. 2d | Fixed; MultiSuSiE capped today |
| 7 | **Disease progression default** (SlopeHunter) never applied because `pipeline.output` always exists; manuscript said SlopeHunter, the demo used CWLS (B27) | Method in text ≠ method run | Code fixed today; **manuscript decision** (report both?) |
| 8 | **Liftover** silently mis-lifted ~0.26% of variants: inverted blocks 2 bp off, minus-strand alleles not complemented, chromosome changes dropped without a report | Small number of wrong positions/alleles after 37→38 | Fixed (liftover QC, drop reasons, report, 5% stop) |
| 9 | **YAML N overrode the per-variant N column** (m19 follow-up) | Wrong N for variants not in every cohort; affects SuSiE scaling and LDSC | Fixed today |
| 10 | **Fine-mapping used the first row's N** for the whole region | As above, per region | Fixed today (median N_EFF/N) |
| 11 | **Z-only input crashed** (`convert_z_score_to_beta` called without N and needing BETA) (B3) | Files with only Z couldn't be used | Fixed today (Zhu et al. 2016 standardised effects) |
| 12 | **I/D indel coding** wasn't put in canonical order, so I/D and D/I gave different keys | The same indel missed between studies | Fixed today |
| 13 | **No strand check** against a reference (B2) | Strand-flipped variants kept with wrong sign | Fixed today (complement non-palindromic, frequency for palindromic, report) |
| 14 | **1000G LD matched on rsID only, no allele check** | r sign wrong at multi-allelic/strand-mismatch sites; LD mismatch in SuSiE | Fixed today (`align_ld_to_gwas`) |
| 15 | **Usage logging** sent user, host, input YAML and error log to the developers' Google sheet whenever the credentials file existed (B35) | Privacy; undocumented | Opt-in today (`GENEHACKMAN_USAGE_LOGGING`) |
| 16 | **Pan-UKB pack LD not positive semi-definite** after quantisation/banding; a diagonal ridge (as first used) distorted SuSiE | New feature, not in submitted results | Fixed: eigenvalue clipping + cov2cor, `LD_MIN_EIGENVALUE` reported |
| 17 | **Clumping r² differs by reference:** PLINK 1.9 uses haplotype r², the Pan-UKB pack genotype-dosage r | Slightly different lead sets between references | Documented; methods sentence in manuscript |
| 18 | **QTL data were GRCh37 while GWAS are lifted to GRCh38**; MetaBrain and trans-QTLs | QTL MR keyed on GRCh37 positions | Replaced: cis-only eQTLGen, BrainMeta, UKB-PPP from genotype-phenotype-map (GRCh38); manuscript text updated |
| 19 | **Fine-mapping window default** was 500 kb in the parser, 1,000 kb in templates (B15) | Inconsistent runs | Unified to 1,000 kb |
| 20 | **Legacy `coloc_analysis()`** treated every trait as quantitative | Not used by the pipelines | Fixed (`coloc_dataset()` handles cc) |

### Not fixed (audit items still open)

| Item | Note |
|---|---|
| B18 | coloc prior sensitivity (report results under other p12) |
| B20 | Steiger filtering; report MR hits that weren't coloc-tested |
| B24–B26 | Winner's curse in expected-vs-observed; Q-test caveats; tidy-ups |
| B28–B31 | Collider-bias method details and labels; `throw()` |
| B34 | e2e tests in CI with output assertions (the known-answer tests are unit tests) |
| B36–B38 | Count skipped loci; DESCRIPTION placeholder; munge `sed` patch |
| C3.4 / D1 | MungeSumstats as an engine (only the concordance script exists) |
| – | Check `subset(mr_results, exposure %in% exposures)` in `R/mr.R` when `exposures` is set (flagged in the handoff, not reproduced) |
| – | Miami/Manhattan titles still use file names (m17) |
| – | LDSC with `nan` prevalences for mixed binary/continuous `--rg` calls: check on ieup1 that LDSC 1.0.1 accepts them (the HF run will show it) |

---

## 3. Draft-letter promises that are not implemented

The tracked response letter now describes what exists. These were removed or reworded; restore them only if they are built before resubmission.

| Draft said | Now says |
|---|---|
| `ldsc.ref_ld_chr` / `ldsc.w_ld_chr` YAML options | `ldsc.ld_scores: pan_ukb` (default) or `1000genomes_b37`; Pan-UKB scores for all ancestries, benchmarked against the Alkes Price scores |
| Alkes Price scores as EUR default | Pan-UKB scores (Andrew's decision, 5 Nov), benchmark reported |
| Precomputed LD matrices (.bcor, .npz) | PLINK genotypes (`ld_bfile`) and the Pan-UKB LD pack |
| Related-work paragraph plus a supplementary table comparing ~14 tools | A few sentences on the different purpose of GeneHackman's standardisation (Andrew, 6 Nov) |
| `finemap.method` selectable | Documented interface and worked example; method option is the extension route |
| MungeSumstats as engine (D1), HyPrColoc (D5), LDAK (D5), Popcorn (D4) | Not added; listed as extensions/limitations |
| Report warns on LDSC-intercept overlap | Static warning in the coloc report and docs; LDSC cross-trait intercept reported |
| TwoSampleMR concordance | MungeSumstats concordance only |
| Identifier parser fills alleles | Fills CHR/BP; alleles come from EA/OA columns |

## 4. What Andrew runs on ieup1 to fill the placeholders

The full step-by-step plan, with commands and the configs in `paper/revision1_hf/`, is `REVISION1_PLAN2.md`. In short:

1. Rebuild the container (EnsDb v86), upload bundle 1.2.0, bump DESCRIPTION.
2. `scripts/validation/compare_with_mungesumstats.R` (see `scripts/validation/README.md`): HFall, BBJ, and BBJ with `--chrpos_only` (Reviewer 1's chr:pos:a1:a2 case). This needs dbSNP and MungeSumstats' reference packages, so it can't run in Claude's sandbox; it is the main outstanding check for R1-M3 and R2-1.
3. LDSC benchmark (runbook step 8) and the HF compare_gwases run with `study_type: binary`, `N_case`, `N_control`, `prevalence` (HF prevalence: decide, e.g. 0.02).
4. HF fine-mapping: `max_causal: auto` (L = 1) and `max_causal: 10`, 1000G and Pan-UKB; collect `LD_MISMATCH_S`, `N_LD_OUTLIERS`, credible-set counts.
5. coloc with N_EFF; Fig. 2d.
6. Disease progression with SlopeHunter (and CWLS for comparison).
7. Clean the Supplementary YAMLs; dataset table.
