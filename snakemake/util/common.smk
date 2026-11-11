import csv
import hashlib
import zipfile
import gzip
import yaml
import os
import re
import sys
import subprocess
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from types import SimpleNamespace

include: "constants.smk"
include: "log_results.smk"


def dict_to_namespace(obj):
    """Recursively convert dict/list from YAML to SimpleNamespace (same shapes as former JSON parsing)."""
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [dict_to_namespace(item) for item in obj]
    return obj


def get_docker_container():
    version = DOCKER_VERSION
    if not version:
        raise ValueError(
            "Error: could not determine Docker image version. "
            "Set DOCKER_VERSION in .env or ensure DESCRIPTION contains Version:"
        )
    pipeline_genomic_dir = os.path.join(PIPELINE_DATA_DIR.rstrip("/"), "genomic_data", "pipeline")
    sif_path = os.path.join(pipeline_genomic_dir, f"genehackman_{version}.sif")
    if not os.path.isfile(sif_path):
        raise ValueError(f"Error: SIF file not found: {sif_path}")
    return sif_path


def parse_pipeline_input(pipeline_includes_clumping=False, allow_flip_false=False, output_builds=None, output_builds_reason=""):
    global input_file
    env_in = os.environ.get("GENEHACKMAN_INPUT") or os.environ.get("GENEHACKMAN_INPUT_FILE")
    try:
        cfg_in = config.get("genehackman_input")
    except NameError:
        cfg_in = None
    if cfg_in is not None and str(cfg_in).startswith("-"):
        cfg_in = None
    input_file = env_in or cfg_in or "input.yaml"

    if not os.path.isfile(".env"):
        raise ValueError("Error: .env file doesn't exist")
    if not os.path.isfile(input_file):
        raise ValueError(f"Error: {input_file} file doesn't exist")

    with open(input_file) as pipeline_input:
        try:
            raw = yaml.safe_load(pipeline_input)
            if raw is None:
                raw = {}
            pipeline = dict_to_namespace(raw)
        except Exception as e:
            raise Exception(
                "ERROR: There is an error with the pipeline YAML file, "
                "please ensure it is valid YAML: https://www.yamllint.com/"
            ) from e

    if not hasattr(pipeline, "is_test"): pipeline.is_test = False
    if not hasattr(pipeline, "output"): pipeline.output = default_output_options
    else:
        if not hasattr(pipeline.output, "build"): pipeline.output.build = default_output_options.build
        elif pipeline.output.build not in build_options:
            raise ValueError(f"Error: {pipeline.output.build} is not one of the valid options: {build_options}")

        if not hasattr(pipeline.output, "columns"): pipeline.output.columns = default_output_options.columns
    if output_builds and pipeline.output.build not in output_builds:
        stop_with_message(f"This pipeline needs output.build: {' or '.join(output_builds)}, not "
                          f"{pipeline.output.build}.\n{output_builds_reason}")
    # ld_reference has no default: the LD source is the user's choice (checked once the GWAS are read, below)
    pipeline.ld_reference = str(getattr(pipeline, "ld_reference", "") or "").strip()
    if pipeline.ld_reference and pipeline.ld_reference not in ld_reference_options:
        stop_with_message(f"ld_reference must be one of {', '.join(ld_reference_options)}, not {pipeline.ld_reference}.")
    if pipeline_includes_clumping:
        clump_arguments = getattr(pipeline, "plink_clump_arguments", None)
        if not clump_arguments:
            raise ValueError(
                'Error: plink_clump_arguments is required for this pipeline, '
                'e.g. plink_clump_arguments: "--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"'
            )
        pipeline.plink_clump_arguments = apply_default_clump_arguments(str(clump_arguments))
    if not hasattr(pipeline, "populate_rsid"): pipeline.populate_rsid = False
    if not hasattr(pipeline, "populate_eaf"): pipeline.populate_eaf = False
    if not hasattr(pipeline, "flip_alleles"): pipeline.flip_alleles = True

    if not hasattr(pipeline, "finemap"):
        pipeline.finemap = SimpleNamespace()
    fm = pipeline.finemap
    if not hasattr(fm, "window_kb"): fm.window_kb = 1000
    # SuSiE L (default 10). It is never changed according to the LD reference: when the GWAS is more than
    # max_n_ratio times larger than the LD reference, the log, finemap_regions.tsv and the reports warn instead
    if not hasattr(fm, "max_causal"): fm.max_causal = 10
    if not str(fm.max_causal).strip().isdigit() or int(str(fm.max_causal).strip()) < 1:
        stop_with_message(f"finemap.max_causal must be a positive whole number (SuSiE L, default 10), not {fm.max_causal}.")
    if not hasattr(fm, "max_n_ratio"): fm.max_n_ratio = 10
    if not hasattr(fm, "coverage"): fm.coverage = 0.95
    if not hasattr(fm, "min_abs_corr"): fm.min_abs_corr = 0.5
    if not hasattr(fm, "max_region_kb"): fm.max_region_kb = 3000
    if not hasattr(fm, "region_strategy"): fm.region_strategy = "merge"
    if fm.region_strategy not in finemap_region_strategies:
        stop_with_message(f"finemap.region_strategy must be one of {', '.join(finemap_region_strategies)}, "
                          f"not {fm.region_strategy}.")
    if not hasattr(fm, "remove_ld_outliers"): fm.remove_ld_outliers = False

    if not hasattr(pipeline, "liftover"):
        pipeline.liftover = SimpleNamespace()
    lo = pipeline.liftover
    if not hasattr(lo, "complement_minus_strand"): lo.complement_minus_strand = False
    if not hasattr(lo, "max_drop_fraction"): lo.max_drop_fraction = 0.05
    if not hasattr(lo, "check_build"): lo.check_build = True

    if not hasattr(pipeline, "harmonise"):
        pipeline.harmonise = SimpleNamespace()
    hm = pipeline.harmonise
    if not hasattr(hm, "strand_check"): hm.strand_check = True
    if not hasattr(hm, "drop_ambiguous_palindromic"): hm.drop_ambiguous_palindromic = False

    if not hasattr(pipeline, "coloc"):
        pipeline.coloc = SimpleNamespace()
    cl = pipeline.coloc
    if not hasattr(cl, "overlap_kb"): cl.overlap_kb = 1000
    if not hasattr(cl, "p1"): cl.p1 = 1e-4
    if not hasattr(cl, "p2"): cl.p2 = 1e-4
    if not hasattr(cl, "p12"): cl.p12 = 5e-6

    for g in pipeline.gwases:
        if not hasattr(g, "N"): g.N = 0
        if not hasattr(g, "columns"): g.columns = SimpleNamespace()
        if not hasattr(g, "remove_extra_columns"): g.remove_extra_columns = False
        if not hasattr(g, "flip_alleles"): g.flip_alleles = getattr(pipeline, "flip_alleles", True)
        if g.file.endswith("_std.tsv.gz"):
            # already standardised: used as is, so it must already be in the output build
            declared = getattr(g, "build", None)
            if declared != pipeline.output.build:
                raise ValueError(
                    f"Error: {g.file} is already standardised, so it is used without liftover. Declare its build "
                    f"(build: {pipeline.output.build}) to confirm it is in the output build {pipeline.output.build}"
                    + (f"; it is declared as {declared}. Re-standardise the original GWAS instead." if declared else ".")
                )
        if not hasattr(g, "build"): g.build = "GRCh37"
        if g.build not in build_options:
            raise ValueError(f"Error: build {g.build} of {g.file} is not one of the valid options: {build_options}")
        resolve_sample_size(g)
        g.label = str(getattr(g, "label", "") or "").strip()
        g.ld_bfile = str(getattr(g, "ld_bfile", "") or "").strip()
        if g.ld_bfile:
            g.ld_bfile = os.path.abspath(g.ld_bfile)
            missing_bfile = [g.ld_bfile + ext for ext in (".bed", ".bim", ".fam") if not os.path.isfile(g.ld_bfile + ext)]
            if missing_bfile:
                stop_with_message(f"ld_bfile of {g.file}: {', '.join(missing_bfile)} not found. ld_bfile is the "
                                  "PLINK prefix of in-sample (or other user-supplied) genotypes, in the output build.")
        if not hasattr(g, "populate_rsid"): g.populate_rsid = False
        if not hasattr(g, "populate_eaf"): g.populate_eaf = False
        g.file = os.path.abspath(g.file)
        g.populate_rsid = resolve_rsid_population(pipeline_includes_clumping, g.populate_rsid or pipeline.populate_rsid)
        g.populate_eaf = bool(g.populate_eaf or pipeline.populate_eaf)
        if hasattr(g, "ancestry") and g.ancestry is not None and str(g.ancestry).strip():
            g.ancestry = str(g.ancestry).strip()
        else:
            g.ancestry = ""
        if g.populate_eaf and not g.ancestry:
            raise ValueError(
                'When populate_eaf is true, each GWAS must include "ancestry" in the input JSON '
                '(same as compare_gwases: AFR, AMR, EAS, EUR, or SAS), e.g. "ancestry": "EUR".'
            )
        g.standardised_memory = estimate_memory_needed_for_standardisation(g.file, g.populate_rsid)
        g.prefix = file_prefix(g.file)
        if not g.label: g.label = g.prefix
        g.vcf_columns = get_columns_for_vcf_parsing(g.columns)
        g.input_columns = resolve_gwas_columns(g.file, g.columns)
        g.output_columns = resolve_gwas_columns(g.file, pipeline.output.columns, check_input_columns=False)
        g.standardised_gwas = standardised_gwas_name(g.file)
        setattr(pipeline,g.prefix,g)
    eaf_ancestries = [g.ancestry for g in pipeline.gwases if g.populate_eaf]
    strand_checked = pipeline.harmonise.strand_check and any(g.ancestry for g in pipeline.gwases)
    if not pipeline.ld_reference:
        if pipeline_includes_clumping or eaf_ancestries or strand_checked:
            stop_with_message(
                "ld_reference is required: choose the LD reference panel for this pipeline.\n"
                "  ld_reference: pan_ukb       precomputed LD from the Pan-UK Biobank (~420,000 EUR; smaller for other "
                "ancestries)\n"
                "  ld_reference: 1000genomes   1000 Genomes high-coverage panels (~500 people per super-population)\n"
                "In-sample genotypes can be given per GWAS with ld_bfile. When a GWAS is much larger than the LD "
                "reference, the pipeline warns that reference LD may miscalibrate fine-mapping.")
        pipeline.ld_reference = "1000genomes"  # not used by this pipeline (no clumping, EAF or strand check)
    if eaf_ancestries:
        validate_ancestries(eaf_ancestries, pipeline.ld_reference)
    if pipeline.ld_reference == "pan_ukb" and (pipeline_includes_clumping or eaf_ancestries):
        ld_ancestries = [g.ancestry for g in pipeline.gwases if pipeline_includes_clumping or g.populate_eaf]
        validate_ancestries(ld_ancestries, pipeline.ld_reference)
    if not allow_flip_false:
        for g in pipeline.gwases:
            if getattr(g, "flip_alleles", True) is False:
                raise ValueError("flip_alleles = FALSE is only available for standardise_gwas.smk")
    return pipeline

def resolve_sample_size(g):
    """Sample size settings of one GWAS.

    study_type: continuous (default) or binary (case-control).  For binary traits N_case and N_control give the
    effective sample size N_eff = 4 / (1/N_case + 1/N_control), and prevalence (the population prevalence) gives
    liability-scale heritability in LDSC.  g.analysis_N is the sample size downstream steps fall back on when the
    standardised GWAS has no per-variant N_EFF or N: N_eff for binary traits, else N."""
    study_type = str(getattr(g, "study_type", "continuous") or "continuous").strip().lower()
    if study_type in ("continuous", "quantitative", "quant"):
        g.study_type = "continuous"
    elif study_type in ("binary", "case_control", "case-control", "categorical", "cc"):
        g.study_type = "binary"
    else:
        stop_with_message(f"study_type of {g.file} must be continuous or binary, not {study_type}.")
    g.N_case = _positive_number(getattr(g, "N_case", getattr(g, "N_cases", None)), "N_case", g.file)
    g.N_control = _positive_number(getattr(g, "N_control", getattr(g, "N_controls", None)), "N_control", g.file)
    g.prevalence = _positive_number(getattr(g, "prevalence", None), "prevalence", g.file)
    if g.prevalence is not None and g.prevalence >= 1:
        stop_with_message(f"prevalence of {g.file} is the population prevalence, a fraction between 0 and 1.")
    g.N_eff = None
    g.sample_prevalence = None
    if g.study_type == "binary":
        if (g.N_case is None) != (g.N_control is None):
            stop_with_message(f"{g.file}: give both N_case and N_control, or neither.")
        if g.N_case is not None:
            g.N_eff = 4 / (1 / g.N_case + 1 / g.N_control)
            g.sample_prevalence = g.N_case / (g.N_case + g.N_control)
            if not getattr(g, "N", 0): g.N = g.N_case + g.N_control
    elif g.N_case is not None or g.N_control is not None or g.prevalence is not None:
        stop_with_message(f"{g.file}: N_case, N_control and prevalence are for binary traits; set study_type: binary.")
    g.analysis_N = g.N_eff if g.N_eff else g.N


def _positive_number(value, name, gwas_file):
    if value is None or str(value).strip() == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        stop_with_message(f"{name} of {gwas_file} must be a number, not {value}.")
    if number <= 0:
        stop_with_message(f"{name} of {gwas_file} must be positive, not {value}.")
    return number


def apply_default_clump_arguments(clump_arguments):
    """Append GeneHackman's default --clump-r2 / --clump-kb when plink_clump_arguments does not set them."""
    provided = clump_arguments.split()
    added = [f"{flag} {value}" for flag, value in default_clump_arguments.items() if flag not in provided]
    if not added:
        return clump_arguments
    print(f"plink_clump_arguments: adding default(s) {' '.join(added)} "
          "(PLINK's own defaults are --clump-r2 0.5 --clump-kb 250)")
    return " ".join([clump_arguments.strip()] + added)


def resolve_rsid_population(pipeline_includes_clumping, populate_rsid):
    if not populate_rsid and pipeline_includes_clumping:
        return populate_rsid_options.partial
    elif not populate_rsid and not pipeline_includes_clumping:
        return populate_rsid_options.none
    else:
        return populate_rsid_options.full


def estimate_memory_needed_for_standardisation(gwas_file, populate_rsid):
    baseline_standardisation_gb = 24
    if populate_rsid != populate_rsid_options.full: return baseline_standardisation_gb

    if gwas_file.endswith(".gz") or gwas_file.endswith(".zip"):
        number_of_lines = subprocess.check_output(f"zcat {gwas_file} | wc -l", shell=True, universal_newlines=True)
    else:
        number_of_lines = subprocess.check_output(f"wc -l < {gwas_file}", shell=True, universal_newlines=True)

    number_of_lines = int(number_of_lines)
    if number_of_lines < 10_000_000:
        return 80
    elif number_of_lines < 35_000_000:
        return 130
    else:
        return 200


def resolve_gwas_columns(gwas_file, column_name_map=None, additional_mandatory_columns=[], check_input_columns=True):
    if isinstance(column_name_map, str):
        column_name_map = read_predefined_column_map(column_name_map)

    if not bool(column_name_map):
        column_name_map = SimpleNamespace()

    column_name_map = vars(column_name_map)
    column_name_map = default_columns | column_name_map

    all_mandatory_columns = list(set(default_mandatory_columns + additional_mandatory_columns))
    mandatory_column_names_in_gwas = [column_name_map[name] for name in all_mandatory_columns]

    ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, column_name_map, check_input_columns)

    cli_string = turn_dict_into_cli_string(column_name_map)
    return cli_string


def read_predefined_column_map(predefined_map_name):
    with open("inst/extdata/predefined_column_maps.csv") as file:
        reader = list(csv.DictReader(file, delimiter=","))
        results = list(filter(lambda x: x['name'] in ("default", predefined_map_name), list(reader)))
        remove_empty = {k: v for k, v in results[1].items() if v}
        return SimpleNamespace(**{**results[0], **remove_empty})


def get_columns_for_vcf_parsing(columns):
    if isinstance(columns, str):
        with open("inst/extdata/predefined_column_maps.csv") as file:
            reader = list(csv.DictReader(file, delimiter=","))
            results = list(filter(lambda x: x['name'] in (columns), list(reader)))
            columns = {k: v for k, v in results[0].items() if v}
            columns = SimpleNamespace(**columns)
        vcf_column_string = ','.join(list(columns.__dict__.values())[1:]),
    else:
        vcf_column_string = ','.join(list(columns.__dict__.values())),
    return vcf_column_string


def ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, column_name_map, check_input_columns):
    if not Path(gwas_file).is_file():
        raise ValueError(f"Error: GWAS file {gwas_file} does not exist")
    if not check_input_columns: return

    if not ".vcf" in gwas_file:
        if gwas_file.endswith(".gz"):
            with gzip.open(gwas_file) as f:
                gwas_headers = f.readline().decode("utf-8").strip()
        elif gwas_file.endswith(".zip"):
            with zipfile.ZipFile(gwas_file) as zf:
                unzipped_file = f"{os.path.dirname(gwas_file)}/{zf.namelist()[0]}"
                with open(unzipped_file) as f:
                    gwas_headers = str(f.readline()).strip()
        else:
            with open(gwas_file) as f:
                gwas_headers = str(f.readline()).strip()

        gwas_headers = re.sub('["\']', '', gwas_headers)
        gwas_headers = re.split('\n|,| |\t',gwas_headers)

        missing = set(mandatory_column_names_in_gwas) - set(gwas_headers)
        if len(missing) > 0:
            raise ValueError(f"Error: {gwas_file} doesn't contain {missing}")

        snp_option_names = [column_name_map.get(name) for name in snp_options]
        missing = set(snp_option_names) - set(gwas_headers)
        if len(missing) == len(snp_option_names):
            raise ValueError(f"Error: {gwas_file} doesn't contain a map for SNP or CHR/BP.  Include one please")

        p_option_names = [column_name_map.get(name) for name in p_options]
        missing = set(p_option_names) - set(gwas_headers)
        if len(missing) == len(p_option_names):
            raise ValueError(f"Error: {gwas_file} doesn't contain a map for P or LOG_P.  Include one please")

        beta_and_or_check = []
        for beta_and_or_option in beta_and_or_options:
            option_column_names = [column_name_map[name] for name in beta_and_or_option]
            missing = set(option_column_names) - set(gwas_headers)
            beta_and_or_check.append(len(missing) > 0)

        if all(beta_and_or_check):
            raise ValueError(f"""Error: {gwas_file} doesn't contain the correct pairings for BETA or OR.
                The options available are: (BETA, SE), or (OR, OR_LB, OR_UB)""")


def stop_with_message(message):
    """Stop with a plain message for the user, without a Python traceback.

    For problems with the input or the downloaded data found while the Snakefile is read, before any job runs
    or any lock is taken (os._exit skips Snakemake's own error report, which would add a traceback)."""
    bar = "=" * 100
    sys.stdout.flush()  # anything printed before goes first
    print(f"\n{bar}\n{message.strip()}\n{bar}\n", file=sys.stderr)
    sys.stderr.flush()
    os._exit(1)


def validate_ancestries(ancestries, ld_reference):
    if ld_reference == "pan_ukb":
        allowed_ancestries = pan_ukb_ancestries + list(pan_ukb_ancestry_aliases)
    else:
        allowed_ancestries = thousand_genomes_ancestries
    unknown = sorted({str(a) for a in ancestries if a not in allowed_ancestries})
    if unknown:
        hint = ""
        if ld_reference != "pan_ukb" and set(unknown) <= set(pan_ukb_ancestries):
            hint = f"\n{', '.join(unknown)} is only available in Pan-UKB: set ld_reference: pan_ukb to use it."
        stop_with_message(f"Ancestry {', '.join(unknown)} can't be used with ld_reference: {ld_reference}.\n"
                          f"Each GWAS's ancestry must be one of: {', '.join(allowed_ancestries)}.{hint}")


def check_reference_data(pipeline, pipeline_includes_clumping, ld_scores=None, ld_score_ancestries=(),
                         qtl_dataset=None):
    """Stop before any job runs if reference data this pipeline needs hasn't been downloaded.

    download_reference_data.sh fetches only the chosen ancestries, dbSNP can be left out, the GRCh37 files
    (output.build: GRCh37) come only with --grch37, and QTL datasets only with --qtl, so a data directory may hold
    only part of the reference data.  Everything missing is listed in one message, with one command that
    downloads it all.

    ld_scores / ld_score_ancestries: the LDSC score set and ancestries (compare_gwases); qtl_dataset: qtl_mr."""
    build = pipeline.output.build
    data_dir = os.path.normpath(PIPELINE_DATA_DIR) if PIPELINE_DATA_DIR else "<PIPELINE_DATA_DIR>"
    used_ancestries = {g.ancestry for g in pipeline.gwases if g.ancestry} | set(ld_score_ancestries)
    missing, missing_qtl, needs_dbsnp, needs_grch37 = [], [], False, False

    if build in THOUSAND_GENOMES_BUILDS:
        grch37 = build == "GRCh37"
        panel_dir = thousand_genomes_dir(build)
        ld_ancestries = {g.ancestry for g in pipeline.gwases if g.ancestry and (pipeline_includes_clumping or g.populate_eaf)}
        if pipeline.ld_reference == "1000genomes":
            for a in sorted(a for a in ld_ancestries if not os.path.isfile(panel_dir + a + ".bed")):
                missing.append(f"the 1000 Genomes {build} LD panel for {a}")
                needs_grch37 |= grch37
        else:
            for a in sorted(ld_ancestries):
                population = pan_ukb_ancestry_aliases.get(a, a)
                if not os.path.isfile(os.path.join(LD_PACK_DIR, population, "info.tsv")):
                    missing.append(f"the Pan-UKB LD pack for {a}" + (f" ({population} in Pan-UKB)" if population != a else ""))
        rsid_options = {g.populate_rsid for g in pipeline.gwases}
        if populate_rsid_options.partial in rsid_options and not os.path.isfile(panel_dir + "marker_to_rsid.tsv.gz"):
            missing.append(f"the 1000 Genomes {build} RSID marker map (to fill in RSIDs)")
            needs_grch37 |= grch37
        dbsnp_dir = GENOMIC_DATA_DIR + "dbsnp/" + THOUSAND_GENOMES_BUILDS[build] + "/"
        if populate_rsid_options.full in rsid_options and not os.path.isfile(dbsnp_dir + "chr1.fst"):
            missing.append(f"dbSNP {build}, used by populate_rsid: true (or leave populate_rsid unset to fill RSIDs "
                           "from the 1000 Genomes marker map instead)")
            needs_dbsnp = True
            needs_grch37 |= grch37
    if ld_scores:
        for a in sorted(set(ld_score_ancestries)):
            population = pan_ukb_ancestry_aliases.get(a, a) if ld_scores == "pan_ukb" else a
            if not os.path.isdir(LDSC_DIRS[ld_scores] + population):
                missing.append(f"the {ld_scores} LD scores for {a} (for LDSC)")
                needs_grch37 |= ld_scores == "1000genomes_b37"
    if qtl_dataset:
        top_hits = (QTL_DIRECTORY or "") + qtl_dataset + "/GRCh38/top_hits/" + qtl_subcategory + ".tsv.gz"
        if not QTL_DIRECTORY or not os.path.isfile(top_hits):
            missing_qtl.append(f"the {qtl_dataset} QTL data")

    if not missing and not missing_qtl:
        return
    nothing_downloaded = not os.path.isdir(GENOMIC_DATA_DIR + "liftover")
    if nothing_downloaded:
        lines = [f"No GeneHackman reference data was found in PIPELINE_DATA_DIR ({data_dir}).",
                 "Download it (only what this pipeline needs) with:"]
    else:
        lines = ["This pipeline needs reference data that hasn't been downloaded yet.", ""]
        if missing:
            lines += [f"Missing from PIPELINE_DATA_DIR ({data_dir}):"] + [f"  - {m}" for m in missing]
        if missing_qtl:
            where = os.path.normpath(QTL_DIRECTORY) if QTL_DIRECTORY else "QTL_DATA_DIR, which isn't set in .env"
            lines += [f"Missing from {('QTL_DATA_DIR (' + where + ')') if QTL_DIRECTORY else where}:"]
            lines += [f"  - {m}" for m in missing_qtl]
        lines += ["", "Download what's missing with:"]

    # one command for everything; the ancestries the pipeline uses (or, if none, those already downloaded)
    ancestries = used_ancestries or {d for d in pan_ukb_ancestries if os.path.isdir(os.path.join(LD_PACK_DIR, d))} or {"EUR"}
    command = f"./download_reference_data.sh --dir {data_dir} --ancestries {','.join(sorted(ancestries))}"
    dbsnp_dir = GENOMIC_DATA_DIR + "dbsnp/" + THOUSAND_GENOMES_BUILDS.get(build, "b38_dbsnp156") + "/"
    if not needs_dbsnp and not os.path.isfile(dbsnp_dir + "chr1.fst"):
        command += " --no-dbsnp"  # don't add dbSNP (11.7 GB) to a directory set up without it
    if needs_grch37:
        command += " --grch37"
    if missing_qtl:
        command += f" --qtl {qtl_dataset}" + (f" --qtl-dir {os.path.normpath(QTL_DIRECTORY)}" if QTL_DIRECTORY else "")
    lines += [f"  {command}", ""]
    lines += ["Then run the pipeline again." if nothing_downloaded else
              "Files you already have are kept, so only the missing ones are downloaded. Then run the pipeline again."]
    if missing_qtl and not QTL_DIRECTORY:
        lines += ["Set QTL_DATA_DIR in .env to the QTL directory the download script prints."]
    lines += ["See ./download_reference_data.sh --help, and 'Download reference data' in README.md."]
    stop_with_message("\n".join(lines))


def standardised_gwas_name(gwas_name):
    if gwas_name.endswith("_std.tsv.gz"):
        return gwas_name
    elif gwas_name.startswith("{") and gwas_name.endswith("}"):
        return DATA_DIR + "gwas/" + gwas_name + "_std.tsv.gz"
    else:
        return DATA_DIR + "gwas/" + file_prefix(gwas_name) + "_std.tsv.gz"

def cleanup_old_slurm_logs():
    if not os.path.isdir(pipeline_log_directory): return

    one_month_ago = datetime.now() - relativedelta(months=1)
    files = [f for f in os.listdir(pipeline_log_directory) if os.path.isfile(os.path.join(pipeline_log_directory, f))]

    for filename in files: 
        file = os.path.join(pipeline_log_directory, filename)
        file_timestamp = datetime.utcfromtimestamp(os.stat(file).st_mtime)
        if file_timestamp < one_month_ago: os.remove(file)




def file_prefix(filename):
    stem = Path(rf"{filename}").stem
    prefix = stem.split('.')[0]
    prefix = re.sub("_std", "", prefix)
    return prefix


def gwas_run_label(gwases):
    """Stable label from GWAS file prefixes for run-specific completion sentinels."""
    labels = sorted({g.prefix for g in gwases})
    label = "_".join(labels)
    if len(label) > 180:
        digest = hashlib.sha1(label.encode("utf-8")).hexdigest()[:12]
        label = f"{len(labels)}gwases_{digest}"
    return label


def multi_finemap_complete_file(gwases):
    label = gwas_run_label(gwases)
    return RESULTS_DIR + f"finemap/multi_ancestry/finemap_complete_{label}.txt"


def locus_zoom_complete_file(gwases):
    label = gwas_run_label(gwases)
    return RESULTS_DIR + f"coloc/locus_zoom/locus_zoom_complete_{label}.txt"


def turn_dict_into_cli_string(results_dict):
    return ','.join(['%s=%s' % (key, value) for (key, value) in results_dict.items()])


def onsuccess(pipeline_name, files_created=list(), results_file=None, is_test=False):
    print("\nPipeline finished, no errors.  List of created files:")
    print(*files_created, sep='\n')

    if results_file:
        print("\n---------------------")
        print("PLEASE VIEW THIS HTML FILE FOR A SUMMARY OF RESULTS:")
        print(f"{results_file}.")
        print(f"\n\n\033[1;35;40m Please do us a favour and cite this pipeline: https://doi.org/10.5281/zenodo.10624713 \033[0m\n")

    update_google_sheet(pipeline_name, succeeded=True, is_test=is_test)


def onerror_message(pipeline_name, is_test=False):
    print("\n---------------------")
    print("There are some documented common errors here: https://github.com/MRCIEU/GeneHackman/wiki")
    log_full_path = None
    if os.path.isdir(pipeline_log_directory):
        try:
            entries = [
                os.path.join(pipeline_log_directory, f)
                for f in os.listdir(pipeline_log_directory)
                if os.path.isfile(os.path.join(pipeline_log_directory, f))
            ]
            if entries:
                log_full_path = max(entries, key=lambda p: os.stat(p).st_mtime)
        except OSError:
            pass
    if log_full_path:
        print("There was an error in the pipeline. If the job ran under Slurm, check the batch log:")
        print(log_full_path)
    else:
        print("There was an error in the pipeline. Check the Snakemake output above,")
        print(f"and logs under: {os.path.join(os.getcwd(), '.snakemake', 'log')}")

    update_google_sheet(pipeline_name, succeeded=False, error_file=log_full_path, is_test=is_test)


cleanup_old_slurm_logs()
