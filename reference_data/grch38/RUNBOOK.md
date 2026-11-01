# Runbook: GRCh38 reference data for GeneHackman 1.2.0 (on ieup1)

The code for GRCh38-first GeneHackman is done; it's on `feature/1.1.1-reviewer-feedback`. This runbook builds the GRCh38 reference data that code expects, checks it, and assembles reference bundle 1.2.0. The plan is in the project doc `claude/grch38_first_plan.md`.

Each step is resumable: re-running a command skips finished work. Run the long steps in tmux.

| Step | What | Time (rough) | Disk |
|---|---|---|---|
| 0 | Second checkout, environment, settings | 30 min | — |
| 1 | Check the dbSNP GRCh38 files | 5 min | — |
| 2 | Build the 1000G 30x GRCh38 panel | 0.5–1 day | ~60 GB temporary, ~10 GB kept |
| 3 | Import the Pan-UKB LD scores | 30 min | ~5 GB |
| 4 | Make GRCh38 copies of the QTL data | hours | as large as the QTL data |
| 5 | Refresh the Pan-UKB pack's GRCh38 keys (after the pack build) | 1–2 h per ancestry | — |
| 6 | Validation report | 30 min | — |
| 7 | Assemble bundle 1.2.0 and the GRCh37 extras | 30 min | hard links, ~0 |
| 8 | LDSC benchmark (reviewer R1-M1) | 1 h | — |
| 9 | End-to-end tests, then release | 1–2 h | — |

When you're done, send me `$WORK_DIR/validation_report.md`, the output of step 3 and the step 8 table.

---

## 0. Set up, without disturbing the running pack build

The Pan-UKB pack build is running from `/home/wt23152/GeneHackman` in tmux. **Don't `git pull` in that checkout while it runs**: bash reads `run_build.sh` as it goes, and changing it under a running script can break the build. Use a second checkout instead:

```bash
cd /home/wt23152/GeneHackman
git fetch origin                     # after pushing feature/1.1.1-reviewer-feedback from the laptop
git worktree add ../GeneHackman_grch38 origin/feature/1.1.1-reviewer-feedback
cd ../GeneHackman_grch38
```

Create the tools environment (plink2, PLINK 1.9, bcftools). R comes from the GeneHackman container, as for the pack build:

```bash
micromamba create -y -p /local-scratch/data/genehackman_grch38/env -f reference_data/grch38/environment.yml
# or: conda env create -p /local-scratch/data/genehackman_grch38/env -f reference_data/grch38/environment.yml
export PATH=/local-scratch/data/genehackman_grch38/env/bin:$PATH
```

Write `reference_data/grch38/local_config.sh`. Git doesn't track it. Fill in the paths you know:

```bash
export WORK_DIR=/local-scratch/data/genehackman_grch38
export OLD_PIPELINE_DATA_DIR=/path/to/genehackman_1.1.0_reference_data   # has genomic_data/ and LDSCORE/
export DBSNP_B38_DIR=/path/to/b38_dbsnp156                              # chr1.fst ... (step 1)
export LD_PACK_DIR=/local-scratch/data/pan_ukb_ld/pack
export OLD_QTL_DATA_DIR=/path/to/genehackman_qtl_1.1.0                  # eqtlgen/, metabrain/
export R_CMD="apptainer exec --bind /local-scratch,/home/wt23152,/path/to/data /local-scratch/data/pan_ukb_ld/genehackman_develop.sif Rscript"
export PATH=/local-scratch/data/genehackman_grch38/env/bin:$PATH
export THREADS=16
```

R loads GeneHackman from this checkout (`GENEHACKMAN_DIR`), so the container only has to provide the R packages. Make sure every directory above is in the apptainer `--bind` list.

**In every new shell**, from `GeneHackman_grch38`, run `source reference_data/grch38/config.sh`. The R and Python steps read their settings from the environment, and apptainer passes it through.

Check that the download links still work. The scripts use these URLs. If one has moved, set the matching variable in `local_config.sh`; the names are in `config.sh`.

```bash
source reference_data/grch38/config.sh
curl -sI "$KG30X_URL/${KG30X_FILE_PATTERN//\{CHR\}/22}" | head -1    # 1000G 30x phased panel, chr22
curl -sI "$KG_PANEL_URL" | head -1                                   # 2,504-person panel file
curl -sI "$PAN_UKB_LDSCORE_URL" | head -1                            # Pan-UKB LD scores
curl -sI https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz | head -1
```

Each should print `HTTP/... 200`.

## 1. Check the dbSNP GRCh38 files

Full RSID population reads `genomic_data/dbsnp/b38_dbsnp156/chr<N>.fst`. The GRCh38 rsIDs for the new 1000G panel come from the same files. They must have the same layout as the GRCh37 ones, which `R/chrpos_to_rsid.R` reads:

```bash
$R_CMD -e '
old <- file.path(Sys.getenv("OLD_PIPELINE_DATA_DIR"), "genomic_data/dbsnp/b37_dbsnp156")
new <- Sys.getenv("DBSNP_B38_DIR")
cat("b37 files:", sort(list.files(old)), "\nb38 files:", sort(list.files(new)), "\n")
for (d in c(old, new)) { f <- fst::fst(file.path(d, "chr22.fst")); print(names(f)); print(head(f[1:3, ])) }'
```

Both should have the same columns (`RSID`, `CHR`, `BP`, `REF`, `ALT`, …) and the same file names, including `chr23`/`chr24` if the b37 set has them. rs4965031 is a quick check: it is at chr22 16052080 in GRCh37 and should sit about 0.5 Mb lower in GRCh38. If the layout differs, stop and send me the output.

## 2. Build the GRCh38 1000 Genomes panel

```bash
tmux new -s grch38
bash reference_data/grch38/10_build_1000g_panel.sh 2>&1 | tee $WORK_DIR/panel.log
```

This builds the panel from the 1000G 30x phased panel (native GRCh38), with the same rules as today's GRCh37 panel:

- the 2,504 unrelated people;
- autosomes;
- biallelic SNVs with MAF ≥ 1% per superpopulation;
- rsIDs from dbSNP 156;
- A1 set to the alphabetical effect allele.

It writes `$WORK_DIR/bundle_1.2.0/genomic_data/1000genomes/b38_dbsnp156/{EUR,EAS,AFR,AMR,SAS}.{bed,bim,fam,frq}` and `marker_to_rsid.tsv.gz`. The log ends with variant and people counts per superpopulation:

- People should be EUR 503, EAS 504, AFR 661, AMR 347 and SAS 489, the same as the GRCh37 panel.
- Variant counts should be close to the GRCh37 panel's, for example EUR about 8.5 M.

Once the panel is checked (step 6), you can delete `$WORK_DIR/1000g_30x/vcf` and `pgen`.

## 3. Import the Pan-UKB LD scores

```bash
source reference_data/grch38/config.sh     # exports the settings the script reads
python reference_data/grch38/20_import_pan_ukb_ldscores.py 2>&1 | tee $WORK_DIR/ldscores.log
```

This writes `bundle_1.2.0/LDSCORE/pan_ukb/<POP>/ld.l2.ldscore.gz` (+ `.M`, `.M_5_50`) for EUR, AFR, AMR, CSA, EAS and MID, plus `w_hm3.snplist`. Pan-UKB doesn't document the archive's layout, so **check the log**:

- it lists what was in the archive;
- for each population, how many LD scores there are;
- whether the SNP column was already rsIDs, or had to be mapped through the LD pack;
- what share of SNPs are in `w_hm3.snplist`. Expect most.

If the script stops because it can't find one file per population, send me the "Archive contents" line.

## 4. GRCh38 copies of the QTL data

```bash
CHAIN_DIR=$OLD_PIPELINE_DATA_DIR/genomic_data/liftover \
  $R_CMD reference_data/grch38/30_lift_qtl_data.R --datasets eqtlgen,metabrain --cores $THREADS
```

For each dataset, the script:

1. works out the build of the files by matching variants against the pack (eQTLGen should be GRCh37; MetaBrain is unknown);
2. lifts GRCh37 files with GeneHackman's liftover and QC, or copies files already in GRCh38;
3. writes `$NEW_QTL_DATA_DIR/<dataset>/GRCh38/{top_hits,gwas}/...`, with a `liftover_summary.tsv`.

Then put the GRCh37 originals beside them as hard links, which use no extra disk:

```bash
for d in eqtlgen metabrain; do cp -al "$OLD_QTL_DATA_DIR/$d" "$NEW_QTL_DATA_DIR/$d/GRCh37"; done
# keep anything else in the old QTL directory (e.g. pqtl/) as it was
for item in "$OLD_QTL_DATA_DIR"/*; do
  name=$(basename "$item"); [[ "$name" == eqtlgen || "$name" == metabrain ]] || cp -al "$item" "$NEW_QTL_DATA_DIR/"
done
```

The new layout is `<dataset>/<build>/{top_hits,gwas}`, which the pipeline now reads using `output.build`.

## 5. Refresh the Pan-UKB pack's GRCh38 keys (after the pack build has finished)

The pack's GRCh38 keys were made with the old liftover, which had no round-trip check. Re-annotating with GeneHackman's new liftover makes them match a GWAS lifted by GeneHackman exactly. The LD files aren't touched. Only once `run_build.sh status` shows every ancestry finished:

```bash
cd /home/wt23152/GeneHackman_grch38          # the new code
bash reference_data/pan_ukb_ld/01_setup_env.sh   # fetches hg38ToHg19 (for the round trip) into the pack's WORK_DIR
bash reference_data/pan_ukb_ld/run_build.sh update-keys     # all ancestries; or: update-keys EUR
bash reference_data/pan_ukb_ld/run_build.sh status
```

`logs/update_keys_<POP>.log` shows the `LIFT_STATUS` counts before and after for each chromosome. Expect `round_trip` to appear, at about 0.03% of variants.

`01_setup_env.sh` re-runs quickly when the environment already exists. If you'd rather not run it, copy `hg38ToHg19.over.chain.gz` into `/local-scratch/data/pan_ukb_ld/liftover_data/genomic_data/liftover/`.

When this has finished, you can switch the original checkout to the new code: `cd /home/wt23152/GeneHackman && git pull`.

## 6. Validation report

```bash
$R_CMD reference_data/grch38/40_validate_grch38.R --pop EUR --chrom 22
$R_CMD reference_data/grch38/40_validate_grch38.R --pop EAS --chrom 22   # optional, other ancestries
```

This writes `$WORK_DIR/validation_report.md`. What to expect:

1. **Pan-UKB GRCh38 keys vs the native GRCh38 panel.** This is the first check of the pack's lifted keys against native GRCh38 data. Most pack keys should be in the panel, with EAF r > 0.95 and almost no |EAF difference| > 0.2.
2. **GeneHackman liftover vs native positions.** The GRCh37 panel is lifted with the new liftover and compared by rsID. Expect close to 100% the same position.
3. **Panel continuity.** r for the same rsID pairs in the old and new panels: correlation > 0.95 and sign agreement close to 100%.
4. **Bundle contents.** Everything present except what step 7 adds.

## 7. Assemble bundle 1.2.0 and the GRCh37 extras

```bash
bash reference_data/grch38/50_assemble_bundle.sh
$R_CMD reference_data/grch38/40_validate_grch38.R    # section 4 should now be all OK
```

`$WORK_DIR/bundle_1.2.0` becomes `gs://genehackman/1.2.0/`. It contains:

- the GRCh38 1000G panel, dbSNP b38, the Pan-UKB pack, the chains including hg18ToHg38, `gene_name_map.tsv`;
- `LDSCORE/pan_ukb`.

`$WORK_DIR/bundle_1.2.0-grch37` becomes `gs://genehackman/1.2.0-grch37/`. It holds the GRCh37 1000G panel, dbSNP b37 and the old LD scores.

**Before uploading**, check the Pan-UKB licence for redistributing the pack and LD scores from the GeneHackman bucket. It is CC BY 4.0, with the note that "usage may be restricted by UK Biobank".

The upload commands, from a machine with write access, are printed at the end. The QTL data go to `gs://genehackman-qtl/1.2.0/`:

```bash
gsutil -m rsync -r $NEW_QTL_DATA_DIR/ gs://genehackman-qtl/1.2.0/
```

## 8. LDSC benchmark (reviewer R1-M1)

Compare h² and rg for the same GWAS with three LD score sets: Pan-UKB (new default), the current GeneHackman 1000G scores (`1000genomes_b37`), and the Alkes Price group's 1000G EUR scores.

1. Download the Alkes Price group's LD scores ("eur_w_ld_chr", per chromosome) and put them at `$WORK_DIR/ldscores_alkes/EUR/`, with `w_hm3.snplist` beside `EUR/`.
2. Use standardised GWAS that have RSIDs and EAF, for example from a `compare_gwases` run on the HF GWAS.

```bash
G1=/path/to/HFall_std.tsv.gz; G2=/path/to/niHFrEF_std.tsv.gz   # standardised GWAS (need RSID, EA, OA, P, BETA, EAF)
NS=N1,N2                                                        # their sample sizes (N_eff for binary traits)
export DATA_DIR=$WORK_DIR/ldsc_benchmark; mkdir -p $DATA_DIR    # run_ldsc.sh writes munged files to $DATA_DIR/ldsc
for set in pan_ukb:$WORK_DIR/bundle_1.2.0/LDSCORE/pan_ukb \
           1000genomes_b37:$OLD_PIPELINE_DATA_DIR/LDSCORE/b37_dbsnp156 \
           alkes:$WORK_DIR/ldscores_alkes; do
  name=${set%%:*}; dir=${set#*:}
  apptainer exec --bind /local-scratch,/path/to/data /local-scratch/data/pan_ukb_ld/genehackman_develop.sif \
    bash scripts/run_ldsc.sh "$G1,$G2" "$NS" EUR "$WORK_DIR/ldsc_benchmark/$name.log" "$dir"
done
grep -A3 "Total Observed scale h2\|Genetic Correlation:" $WORK_DIR/ldsc_benchmark/*.log
```

`run_ldsc.sh` reads either per-chromosome LD scores (`EUR/1.l2.ldscore.gz` …) or the genome-wide file `EUR/ld.l2.ldscore.gz`. Send me the h², intercept and rg for each set; I'll turn them into the supplementary benchmark table.

## 9. End-to-end tests, then release

The end-to-end tests run every pipeline with both LD references, in GRCh38 (default) and GRCh37. So the test data directory needs the main bundle, the GRCh37 extras and the GRCh38 QTL data. Point `.env` at them:

```bash
mkdir -p $WORK_DIR/test_pipeline_data
cp -al $WORK_DIR/bundle_1.2.0/. $WORK_DIR/test_pipeline_data/
cp -aln $WORK_DIR/bundle_1.2.0-grch37/. $WORK_DIR/test_pipeline_data/
# .env: PIPELINE_DATA_DIR=$WORK_DIR/test_pipeline_data  QTL_DATA_DIR=$NEW_QTL_DATA_DIR  DOCKER_VERSION=<image with EnsDb v86>
./tests/e2e_tests/run_test_pipelines.sh
```

The container needs `EnsDb.Hsapiens.v86`, the GRCh38 gene annotation for LocusZoom (added to `docker/requirements.R`). Rebuild the image before running coloc end to end:

```bash
docker build -t mrcieu/genehackman:1.2.0 -f docker/Dockerfile .
```

Then, for the release:

1. Bump `Version:` in `DESCRIPTION` to 1.2.0. The image tag and bucket prefix follow it.
2. Push the image.
3. Upload the bundles (step 7).
4. Commit `tests/testing_complete.txt`.

## If something goes wrong

- **Any step:** re-run the same command; finished work is skipped. To redo a 1000G step, delete its marker in `$WORK_DIR/1000g_30x/state/`.
- **`update-keys` fails** with "pack rows and annotated variants differ": the annotation no longer lines up with the pack. Don't edit anything; send me the log.
- **Standardising a GWAS stops** with "declared as GRCh37, but its positions match GRCh38": the GWAS's `build:` is wrong. That's the new check working.
- **Liftover stops** with "dropped more than 5%": the input build is probably wrong. If it really is right, raise `liftover.max_drop_fraction` in the YAML.
