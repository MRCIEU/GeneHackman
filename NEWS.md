# GeneHackman 1.2.0

Reference data release 1.2.0 (GRCh38), container `mrcieu/genehackman:1.2.0`. This release answers the reviewers'
comments on the manuscript and corrects the errors found while validating the pipelines. Existing YAML files may
need changes: see **Changes that need action**.

## Changes that need action

- **`ld_reference` is required** (`pan_ukb` or `1000genomes`) in every pipeline that clumps, fine-maps, checks
  strand or fills allele frequencies. It has no default, so the choice of LD reference is always explicit.
- **`finemap.max_causal: auto` is removed.** The number of signals per region (SuSiE L) defaults to 10 and is never
  changed by the pipeline. When the GWAS is more than `finemap.max_n_ratio` (default 10) times larger than the LD
  reference, a warning is written to the log, to `finemap_regions.tsv` (`LD_SIZE_WARNING`) and to the HTML reports.
- **GRCh38 is the default output build** (`output.build`). GRCh37 output needs the GRCh37 files
  (`download_reference_data.sh --grch37`). An input file that is already standardised (`*_std.tsv.gz`) must declare
  its build, and it must match the output build.
- **Reference data 1.2.0.** Download it with `./download_reference_data.sh`; release 1.1.0 stays available for
  reproducing earlier results.
- **Usage logging is off** unless `GENEHACKMAN_USAGE_LOGGING=true` is set in `.env`.

## New

- **LD for fine-mapping:** in-sample or other user-supplied genotypes per GWAS (`ld_bfile`), used by SuSiE and
  MultiSuSiE; precomputed Pan-UK Biobank LD (`ld_reference: pan_ukb`).
- **LD diagnostics:** each region is checked for inconsistency between z-scores and LD (susieR `estimate_s_rss`
  and `kriging_rss`; a Python port for MultiSuSiE). Credible sets carry an `LD_FLAG`, which the colocalization and
  QTL reports also show. LD is oriented to the GWAS effect allele, and variants with mismatched alleles are removed.
- **Fine-mapping output:** an HTML report, per-variant PIPs, and region and credible-set tables; regions by
  merging (capped at `max_region_kb`), no merging, or tiling (`finemap.region_strategy`), the same for single- and
  multi-ancestry fine-mapping; regions overlapping the MHC are flagged.
- **LDSC:** Pan-UK Biobank LD scores by default (`ldsc.ld_scores: pan_ukb`); every pair of GWAS of the same
  ancestry; liability-scale heritability for binary traits; genetic correlations flagged when a heritability
  z-score is below 4.
- **Binary traits and sample size:** `study_type`, `N_case`, `N_control` and `prevalence` per GWAS; effective
  sample size, per variant where the file has case and total counts; per-variant N takes precedence over the YAML.
- **Harmonisation:** strand and palindromic checks against the LD reference panel, with a per-GWAS report;
  build conversion that removes and reports unreliable variants (chromosome change, alternative contig, minus strand,
  failed round trip) and stops if more than 5% are removed; a check that the declared input build matches the
  positions; GRCh36 (hg18) input; rsIDs from the dbSNP build of the output build; chr:pos:a1:a2 identifiers; z-score-
  only input; a GWAS Catalog preset.
- **QTL data:** cis-QTLs from eQTLGen, BrainMeta and the UK Biobank Pharma Proteomics Project, in GRCh38, from the
  genotype-phenotype map.
- **Labels:** each GWAS can have a `label`, used in plots and reports.
- **Reproducibility:** each HTML report records the container image, the versions of R and the main packages, and
  the reference-data release.
- **Tests and validation:** known-answer tests for each harmonization rule, simulated fine-mapping tests, and
  `scripts/validation/compare_with_mungesumstats.R`, which compares GeneHackman's output with MungeSumstats.
- **Documentation:** how to add a method to an existing pipeline (CONTRIBUTING.md), and how to run GeneHackman in
  a trusted research environment (README.md).

## Fixed

- The OpenGWAS VCF preset assigned the effect to the reference allele instead of the alternative allele.
- LD scores in reference data 1.1.0 had been computed with a 1-kb window instead of 1 cM, which inflated heritability
  estimates and intercepts. They are no longer used by default.
- LDSC estimated only the genetic correlations of the first GWAS with each of the others.
- Expected-versus-observed replication counted observed replications with a one-sided P value against a two-sided
  expectation, and with four or more GWAS some pairwise results were overwritten.
- Colocalization gave variants without fine-mapping results a log Bayes factor of zero instead of excluding them.
- QTL colocalization used the GWAS sample size and trait type instead of those of the QTL data, and ignored
  user-specified priors.
- When only a P-value threshold was given for clumping, PLINK's defaults (r² 0.5, 250 kb) applied; GeneHackman now
  adds `--clump-r2 0.01 --clump-kb 1000`. PLINK errors now stop the pipeline instead of giving empty results.
- The default index-event bias correction (SlopeHunter) was not applied when no method was specified.
- Build conversion put positions in inverted segments 2 bp off, did not complement minus-strand alleles, and dropped
  variants that moved to another chromosome without a report.
- A sample size given in the YAML overrode the per-variant N in the file, and fine-mapping used the first variant's
  N for a whole region.
- Insertions and deletions coded I/D got a different identifier depending on which allele was the effect allele.
- LD from 1000 Genomes was matched on rsID without checking alleles.
- Multi-ancestry fine-mapping merged neighbouring windows without a size limit.
- GRCh38 files without rsIDs could not be standardised when variants were identified by chr:pos.
