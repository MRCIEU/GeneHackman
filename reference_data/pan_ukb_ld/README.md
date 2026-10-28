# Building the Pan-UKB LD pack

These scripts turn the Pan-UK Biobank LD matrices (precomputed r, no individual-level data, 43 TB, GRCh37, ±10 Mb) into the compact GeneHackman LD pack. The pack keeps:

- variants with MAF ≥ 1% in each ancestry;
- r for variant pairs within ±500 kb;
- r oriented to GeneHackman's effect allele, with keys in both GRCh37 and GRCh38;
- r stored as int8 (r × 127), with |r| < 0.02 stored as 0, xz-compressed.

`PACK_README.md` describes the finished pack. GeneHackman reads it through `R/ld_reference.R` and `scripts/ld_reference.py` when a pipeline sets `ld_reference: pan_ukb`.

The build is meant to run on **ieup1**, in `WORK_DIR=/local-scratch/data/pan_ukb_ld` (the default). It runs on any Linux or macOS machine with enough disk if you set `WORK_DIR`.

## How it saves disk and time

- **Only the blocks that matter are downloaded.** For each ancestry, `05_fetch_source.py` reads the Pan-UKB block index and downloads only the blocks holding r for kept variants within the radius. It then writes a local index listing just those blocks, so Hail reads them as a sparse matrix. The download is roughly the near-diagonal tenth of each ancestry's matrix; the phase 0 report gives the measured figure (`est. download TB`).
- **One ancestry at a time.** An ancestry's blocks are deleted once its pack is written and validated, before the next ancestry starts. Set `KEEP_SOURCE_BLOCKS=1` to keep them. The download refuses to start if it would leave less than `MIN_FREE_GB` free.
- **Everything resumes.** Downloads skip complete files. LD extraction checkpoints every chunk of 2,048 rows. Finished stages leave markers in `state/`. After a crash, reboot or `Ctrl-C`, run the same command again.

## Steps on ieup1

1. **Get the scripts onto ieup1.** Check out the branch with this directory in the GeneHackman clone there, or copy `reference_data/pan_ukb_ld/` and `scripts/ld_reference.py` into it. The R stage loads GeneHackman from that clone (`GENEHACKMAN_DIR`, by default the clone the scripts are in).

2. **Check the machine:**

   ```bash
   bash reference_data/pan_ukb_ld/00_check_environment.sh
   ```

   Look for: free space, tmux, conda/apptainer/R, and an HTTP `200` from the Pan-UKB bucket.

3. **Set up**, once:

   ```bash
   bash reference_data/pan_ukb_ld/01_setup_env.sh
   ```

   This creates `$WORK_DIR/env` (Python 3.11, Java 11, Hail 0.2.139, boto3, UCSC liftOver) and fetches the hg19→hg38 chain. The chain is copied from `$PIPELINE_DATA_DIR/genomic_data/liftover/` when that is set, so the pack uses exactly GeneHackman's chain.

   If plain `Rscript` can't load GeneHackman, point `R_CMD` at the pipeline container in `reference_data/pan_ukb_ld/local_config.sh`:

   ```bash
   export R_CMD="apptainer exec --bind /local-scratch,$HOME /path/to/genehackman_1.1.0.sif Rscript"
   export PIPELINE_DATA_DIR=/path/to/pipeline_data   # optional: reuse GeneHackman's chain file
   ```

4. **Phase 0** (hours, not days):

   ```bash
   reference_data/pan_ukb_ld/run_in_tmux.sh phase0
   ```

   For every ancestry this:
   - downloads the variant and LD-score tables;
   - lists the kept variants and gives them GRCh37 and GRCh38 keys;
   - builds chromosome 22 in full and checks it byte for byte against the source;
   - measures compressed sizes under alternative encodings.

   For EUR it also extracts the test region `6:35.6–37.6 Mb` at ±2 Mb to measure LD beyond 500 kb. If you set `PHASE0_GWAS` to a GeneHackman-standardised GRCh37 GWAS with a signal there, it also compares SuSiE results with full and truncated LD.

   The report, `phase0/report.md`, gives per ancestry: individuals, kept variants, AF coverage, liftover outcomes, estimated pack size, download size and extraction time.

5. **Check alleles against 1000G** (anywhere with PLINK and the 1000G files, e.g. after copying `pack/` to the laptop):

   ```bash
   python reference_data/pan_ukb_ld/61_compare_1000g.py --pack <pack> --pop EUR --chrom 22 \
     --bfile <genomic_data>/1000genomes/b37_dbsnp156/EUR --frq <...>/EUR.frq --window 20000000-25000000
   ```

   Expect EAF correlation ≈ 1 with no allele swaps, and r sign concordance ≥ 98% for |r| > 0.3.

6. **Full build** (days; one ancestry at a time; EUR first gives the most useful pack soonest):

   ```bash
   reference_data/pan_ukb_ld/run_in_tmux.sh full EUR
   reference_data/pan_ukb_ld/run_in_tmux.sh full          # the rest of $POPS
   ```

   Chromosome 22 from phase 0 is reused.

7. **Copy the pack to the laptop**, run from the laptop:

   ```bash
   rsync -av --progress ieup1:/local-scratch/data/pan_ukb_ld/pack/ \
     /Users/wt23152/Documents/Projects/scratch/020/data/genomic_data/pan_ukb_ld/
   ```

## Watching and resuming

| Task | Command |
|---|---|
| Attach to the run | `tmux attach -t ldpack` (detach with `Ctrl-b d`; window `status` refreshes every minute) |
| Progress, disk and last log lines | `reference_data/pan_ukb_ld/run_build.sh status` |
| Logs | `$WORK_DIR/logs/`: `run_build_<date>.log`, per-stage `<POP>_<stage>.log`, Hail logs |
| Resume | start the same `run_in_tmux.sh …` command again (if the old tmux session is still there, `tmux kill-session -t ldpack` first) |

**Changing settings after a stage has run** (for example `RADIUS_KB`): the stage stops with a message naming the marker and outputs to remove. It never silently mixes settings.

## Settings (config.sh)

Every setting can be overridden with an environment variable or in `local_config.sh`, which git doesn't track.

| Setting | Default | Meaning |
|---|---|---|
| `WORK_DIR` | `/local-scratch/data/pan_ukb_ld` | everything is written here |
| `POPS` | `EUR AFR CSA EAS MID AMR` | ancestries to build |
| `MAF_MIN` | 0.01 | minimum MAF |
| `RADIUS_KB` | 500 | radius of stored pairs |
| `DEADZONE` | 0.02 | \|r\| below this is stored as 0 |
| `CHUNK_ROWS` | 2048 | rows per compressed chunk |
| `FETCH_ROWS` | 8192 | rows per source read |
| `HAIL_CORES` | 8 | Hail cores |
| `HAIL_DRIVER_MEMORY` | 32g | Hail driver memory |
| `DOWNLOAD_THREADS` | 16 | parallel downloads |
| `MIN_FREE_GB` | 500 | free space to leave when downloading |
| `KEEP_SOURCE_BLOCKS` | 0 | 1 keeps downloaded blocks after an ancestry is done |

## Files

| File | Stage |
|---|---|
| `05_fetch_source.py` | download tables; plan or download the needed LD blocks; delete blocks |
| `10_extract_variants.py` | Hail: join the variant index to AF, keep biallelic autosomal MAF ≥ 1% variants |
| `20_annotate_variants.R` | GeneHackman `standardise_alleles()` for EA/OA/EAF/SNP37; liftOver to GRCh38 (strand-aware, cross-checked against `convert_reference_build_via_liftover()`) for SNP38 |
| `30_extract_ld.py` | Hail: read r, cut the ±radius band, orient, quantise, compress, checkpoint; `dense` for the phase 0 region |
| `40_finalise_pack.py` | `info.tsv`, `MANIFEST.json` with checksums, the pack README |
| `50_validate_pack.py` | structure checks, and sampled chunks byte-identical to re-reading the source |
| `60_phase0_report.py`, `61_compare_1000g.py`, `62_truncation_finemap.R` | phase 0 measurements |
| `ldpack_writer.py` | the pack writer (the reader is `scripts/ld_reference.py`) |
| `make_test_pack.py` | builds the tiny fixture in `tests/testthat/data/ld_pack_tiny` |
| `tests/run_synthetic_build.sh` | end-to-end test on a synthetic Pan-UKB-format source, including an interrupted and resumed run |
