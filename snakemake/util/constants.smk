from datetime import datetime
from dotenv import load_dotenv
from types import SimpleNamespace
from enum import Enum
import os


def format_dir_string(directory):
    if not directory: return None
    return directory + "/" if not directory.endswith('/') else directory

load_dotenv()
docker_repo = "docker://mrcieu/genehackman"
user = os.getenv('USER')
input_file = "input.yaml"
start_time = datetime.now()

# Slurm writes batch stdout to this directory on HPC; local/Docker runs use ./pipeline_logs instead.
_pipeline_log_override = os.getenv("PIPELINE_LOG_DIR")
_local_flag = os.getenv("GENEHACKMAN_LOCAL", "").strip().lower() in ("1", "true", "yes")
_hpc_work_base = f"/user/work/{user}" if user else None
if _pipeline_log_override:
    pipeline_log_directory = format_dir_string(_pipeline_log_override)
elif _local_flag or not (_hpc_work_base and os.path.isdir(_hpc_work_base)):
    pipeline_log_directory = format_dir_string(os.path.join(os.getcwd(), "pipeline_logs"))
else:
    pipeline_log_directory = format_dir_string(f"{_hpc_work_base}/slurm_logs")

slurm_log_directory = pipeline_log_directory

default_clump_headers = "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2"
#TODO: this should be read from predefined_column_map.csv
default_columns = dict(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA",
    SE="SE", OR="OR", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID", N="N",
    ENSEMBL_ID="ENSEMBL_ID", GENE_NAME="GENE_NAME"
)
#populate_rsid_options = dict(none="none", partial="partial", full="full")
class populate_rsid_options(Enum):
    none = "none"
    partial = "partial"
    full = "full"

default_mandatory_columns = ["EA", "OA"]
snp_options = ["SNP", "BP"]
effect_options = ["BETA", "OR", "Z"]
build_options = ["GRCh38", "GRCh37", "GRCh36"]

p_options = ["P", "LOG_P"]
beta_and_or_options = [
    ["BETA", "SE"],
    ["OR", "OR_LB", "OR_UB"]
]
default_output_options = SimpleNamespace(**{
    "build": "GRCh37",
    "columns": SimpleNamespace(**default_columns)
})

if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")
if not os.getenv("RDFS_DIR"):
    print("Please populate RDFS_DIR in .env if you want the generated files to be automatically copied to RDFS")

DOCKER_VERSION = os.getenv('DOCKER_VERSION')
DATA_DIR = format_dir_string(os.getenv('DATA_DIR'))
RESULTS_DIR = format_dir_string(os.getenv('RESULTS_DIR'))
RDFS_DIR = format_dir_string(os.getenv('RDFS_DIR'))

PIPELINE_DATA_DIR = format_dir_string(os.getenv('PIPELINE_DATA_DIR'))
QTL_DATA_DIR = os.getenv("QTL_DATA_DIR", "").strip()

LDSC_DIR = format_dir_string(PIPELINE_DATA_DIR + "/LDSCORE/b37_dbsnp156")
GENOMIC_DATA_DIR = format_dir_string(PIPELINE_DATA_DIR + "/genomic_data")
THOUSAND_GENOMES_DIR = format_dir_string(PIPELINE_DATA_DIR + "/1000genomes/b37_dbsnp156")
QTL_DIRECTORY = format_dir_string(QTL_DATA_DIR)

if RDFS_DIR and not RDFS_DIR.endswith("working/"):
    raise ValueError("Please ensure RDFS_DIR ends with working/ to ensure the data gets copied to the correct place")
if DATA_DIR == RESULTS_DIR:
    raise ValueError("DATA_DIR and RESULTS_DIR must be different directories")
