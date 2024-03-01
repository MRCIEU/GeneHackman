import csv
import gzip
import json
import os
import time
import re
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from types import SimpleNamespace

include: "constants.smk"
include: "log_results.smk"

def get_docker_container():
    version = "latest"
    with open("../../DESCRIPTION") as file:
        for line in file:
            match = re.match(r"^Version: (\w+)", line)
            if match:
                version = match
                break
    return docker_repo + ":" + version

def parse_pipeline_input():
    if not os.path.isfile(".env"):
        raise ValueError("Error: .env file doesn't exist")
    if not os.path.isfile(input_file):
        raise ValueError(f"Error: {input_file} file doesn't exist")

    with open(input_file) as pipeline_input:
        pipeline = json.load(pipeline_input,object_hook=lambda data: SimpleNamespace(**data))

    if not hasattr(pipeline, "is_test"): pipeline.is_test = False
    if not hasattr(pipeline, "output"): pipeline.output = default_output_options
    else:
        if not hasattr(pipeline.output, "build"): pipeline.output.build = default_output_options.build
        elif pipeline.output.build not in build_options:
            raise ValueError(f"Error: {pipeline.output.build} is not one of the valid options: {build_options}")

        if not hasattr(pipeline.output, "columns"): pipeline.output.columns = default_output_options.columns
    if not hasattr(pipeline, "populate_rsid"): pipeline.populate_rsid = False

    for g in pipeline.gwases:
        if not hasattr(g, "N"): g.N = 0
        if not hasattr(g, "build"): g.build = "GRCh37"
        if not hasattr(g, "populate_rsid"): g.populate_rsid = False
        g.prefix = file_prefix(g.file)
        g.vcf_columns = get_columns_for_vcf_parsing(g.columns)
        g.input_columns = resolve_gwas_columns(g.file, g.columns)
        g.output_columns = resolve_gwas_columns(g.file, pipeline.output.columns, check_input_columns=False)
        g.standardised_gwas = standardised_gwas_name(g.file)
        g.standardised_memory = 56*(g.populate_rsid or pipeline.populate_rsid) + 16
        setattr(pipeline,g.prefix,g)
    return pipeline


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
    return vcf_column_string


def ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, column_name_map, check_input_columns):
    if not Path(gwas_file).is_file():
        raise ValueError(f"Error: {gwas_file} does not exist")
    if not check_input_columns: return

    if not ".vcf" in gwas_file:
        if gwas_file.endswith(".gz"):
            with gzip.open(gwas_file) as f:
                gwas_headers = f.readline().decode("utf-8").strip()
        else:
            with open(gwas_file) as f:
                gwas_headers = str(f.readline()).strip()
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


def validate_ancestries(ancestries):
    allowed_ancestries = ["EUR", "EAS", "AFR", "AMR", "SAS"]
    if not all(ancestry in allowed_ancestries for ancestry in ancestries):
        raise ValueError(f"Please ensure all ancestries are one of these values:\n {allowed_ancestries}")


def standardised_gwas_name(gwas_name):
    if gwas_name.endswith("_std.tsv.gz"):
        return gwas_name
    elif gwas_name.startswith("{") and gwas_name.endswith("}"):
        return DATA_DIR + "gwas/" + gwas_name + "_std.tsv.gz"
    else:
        return DATA_DIR + "gwas/" + file_prefix(gwas_name) + "_std.tsv.gz"

def cleanup_old_slurm_logs():
    if not os.path.isdir(slurm_log_directory): return

    one_month_ago = datetime.now() - relativedelta(months=1)
    files = [f for f in os.listdir(slurm_log_directory) if os.path.isfile(f)]
    print("deleting old logs")

    for filename in files: 
        file = os.path.join(slurm_log_directory, filename)
        file_timestamp = datetime.utcfromtimestamp(os.stat(file).st_mtime)
        if file_timestamp < one_month_ago: os.remove(file)


def file_prefix(filename):
    stem = Path(rf"{filename}").stem
    prefix = stem.split('.')[0]
    prefix = re.sub("_std", "", prefix)
    return prefix


def turn_dict_into_cli_string(results_dict):
    return ','.join(['%s=%s' % (key, value) for (key, value) in results_dict.items()])


def copy_data_to_rdfs(files_created):
    if RDFS_DIR is not None:
        for file_created in files_created:
            rdfs_file = None
            if file_created.startswith(DATA_DIR):
                rdfs_file = file_created.replace(DATA_DIR, RDFS_DIR + "/data/")
            elif file_created.startswith(RESULTS_DIR):
                rdfs_file = file_created.replace(RESULTS_DIR, RDFS_DIR + "/results/")

            if rdfs_file:
                os.system(f"install -D {file_created} {rdfs_file}")
        print(f"Files successfully copied to {RDFS_DIR}")


def onsuccess(pipeline_name, files_created=list(), results_file=None, is_test=False):
    print("\nPipeline finished, no errors.  List of created files:")
    print(*files_created, sep='\n')

    if results_file:
        print("\n---------------------")
        print("PLEASE VIEW THIS HTML FILE FOR A SUMMARY OF RESULTS:")
        print(f"scp {user}@bc4login1.acrc.bris.ac.uk:{results_file} .")

    copy_data_to_rdfs(files_created)

    update_google_sheet(pipeline_name, succeeded=True, is_test=is_test)


def onerror_message(pipeline_name, is_test=False):
    last_log = subprocess.check_output(f"ls -t {slurm_log_directory} | head -n1", shell=True, universal_newlines=True)
    log_full_path = slurm_log_directory + last_log
    print("\n---------------------")
    print("There was an error in the pipeline, please check the last written slurm log to see the error:")
    print(log_full_path)

    update_google_sheet(pipeline_name, succeeded=False, error_file=log_full_path, is_test=is_test)


cleanup_old_slurm_logs()
