# GeneHackman LD pack: Pan-UKB, MAF ≥ 1%, ±1 Mb

Precomputed linkage disequilibrium (Pearson r between variants) for six ancestry groups. The pack holds no individual-level data. GeneHackman reads it when a pipeline sets `ld_reference: pan_ukb`.

## Source

The LD comes from the Pan-UK Biobank LD matrices (`UKBB.<POP>.ldadj.bm`): <https://pan.ukbb.broadinstitute.org/docs/ld>.

- **Licence:** CC BY 4.0, with the note that "usage may be restricted by UK Biobank".
- **Attribution:** please cite the Pan-UKB team (Karczewski et al., *Nature Genetics* 2025) when you use these data.
- **Build provenance:** `MANIFEST.json` records the build settings, the GeneHackman commit and sha256 checksums.

## Contents

Per ancestry, in `<POP>/`:

| File | Contents |
|---|---|
| `info.tsv` | format version, int8 scale, radius, MAF threshold, dead-zone, number of individuals the LD was estimated from |
| `variants_chr<C>.tsv.gz` | one row per variant, in GRCh37 order (columns below) |
| `ld_chr<C>.bin`, `ld_chr<C>.chunks.tsv` | r for every variant pair within the radius, as int8 (r × 127), in xz-compressed chunks of rows |

Columns of `variants_chr<C>.tsv.gz`:

- `IDX`: row index.
- `CHR`, `BP37`, `BP38`: chromosome and positions in each build.
- `EA`, `OA`: effect and other allele, with EA the alphabetically smaller allele, as in GeneHackman's standardised GWAS.
- `EAF`: frequency of EA.
- `SNP37`, `SNP38`: the GeneHackman keys `CHR:BP_EA_OA`.
- `SEG_LEN`: the number of later variants within the radius.
- `LIFT_STATUS`: `ok`, `unmapped`, `strand_flip`, `chrom_change`, `duplicate_grch38` or `duplicate_grch37`. When it isn't `ok`, `BP38` and `SNP38` are empty. `duplicate_grch37` marks variants that share a GeneHackman key with another variant (e.g. an insertion and a deletion with the same alleles at one position); they have no `SNP37` either, because a GWAS row with that key could be either variant.
- `REF37`, `ALT37`, `RSID_UKB`: the source alleles and the rsID given by Pan-UKB.

r is oriented to EA, so it can be used directly with the z-scores of a GeneHackman-standardised GWAS in either build. Pairs further apart than the radius are not stored, and |r| below the dead-zone is stored as 0.

The pack ancestries are EUR, AFR, AMR, CSA, EAS and MID. GeneHackman maps the ancestry code SAS to CSA.

Reader code: `R/ld_reference.R` and `scripts/ld_reference.py` in the GeneHackman repository. How the pack is built: `reference_data/pan_ukb_ld/README.md`.
