# disease_progression

Compare an **incident** and **subsequent** GWAS for collider bias, apply a correction method, and report Miami plots plus expected-vs-observed metrics before and after adjustment.

**Workflow file:** [`disease_progression.smk`](disease_progression.smk)  
**Example input:** [`input_templates/disease_progression.yaml`](input_templates/disease_progression.yaml)

## Run

```bash
./run_pipeline.sh snakemake/disease_progression.smk path/to/input.yaml
```

## Input

### GWAS count

- **Exactly two** GWAS under `gwases:`:
  1. **Incident** GWAS (first entry)
  2. **Subsequent** GWAS (second entry)

### Per-GWAS fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | Yes | Input summary-statistics file. |
| `ancestry` | **Yes** | Used for clumping LD reference. May differ between incident and subsequent. |
| `columns` | No | Column map or preset. |
| `build` | No | Default `GRCh37`. |
| `populate_rsid` | No | See [PIPELINES.md](../PIPELINES.md#shared-yaml-schema). |

### Root-level fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | PLINK clump settings (clumping uses the **incident** GWAS leads). |
| `output.adjusted_gwas.type` | No | Correction method: `slopehunter`, `cwls`, or `mr_ivw`. Default `slopehunter`. |
| `output.adjusted_gwas.p_val` | No | P-value threshold for the adjustment. Default `0.001`. Allowed: `0.1`, `0.01`, `0.001`, `1e-5`. |
| `populate_rsid` | No | Default `false`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise both GWAS
2. Clump incident and subsequent GWAS separately
3. Collider-bias correction on subsequent GWAS (`correct_for_collider_bias.R`) using incident clumps
4. Miami plots: incident vs adjusted; subsequent vs adjusted
5. Expected vs observed for (incident, subsequent) and (incident, adjusted)
6. HTML report

## Output

### Data (`data/`)

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | Standardised incident and subsequent GWAS. |
| `clumped_snps/<prefix>.clumped` | Clumped files for both GWAS. |

### Results (`results/`)

| Path | Description |
| ---- | ----------- |
| `collider_bias/<subsequent_prefix>_collider_bias_results.tsv` | Collider bias analysis results. |
| `collider_bias/<subsequent_prefix>_harmonised_effects.tsv.gz` | Harmonised effect estimates. |
| `collider_bias/<subsequent_prefix>_<type>_adjusted.tsv.gz` | Collider-adjusted subsequent GWAS (e.g. `_cwls_adjusted.tsv.gz`). |
| `collider_bias/original_expected_vs_observed_outcomes.tsv` | Expected vs observed (incident vs subsequent). |
| `collider_bias/original_expected_vs_observed_variants.tsv` | Variant-level replication (original pair). |
| `collider_bias/adjusted_expected_vs_observed_outcomes.tsv` | Expected vs observed (incident vs adjusted). |
| `collider_bias/adjusted_expected_vs_observed_variants.tsv` | Variant-level replication (adjusted pair). |
| `collider_bias/result_<incident_prefix>_<subsequent_prefix>.html` | HTML report. |
| `plots/<subsequent_prefix>_miami_plot.png` | Miami plot (incident vs adjusted). |
| `plots/<adjusted_prefix>_miami_plot.png` | Miami plot (subsequent vs adjusted). |

## Notes

- The HTML report may include guidance on refitting the subsequent GWAS from collider-bias output.
- Clumping for expected-vs-observed comparisons uses the **incident** clumped file only.
