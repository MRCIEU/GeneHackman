from datetime import datetime
from dotenv import load_dotenv
from pathlib import Path
from types import SimpleNamespace
from enum import Enum
import os


def format_dir_string(directory):
    if not directory: return None
    return directory + "/" if not directory.endswith('/') else directory


def _available_cpus():
    """CPUs on the allocated Slurm node if Slurm sets them, else logical CPUs on this host."""
    slurm = os.getenv("SLURM_CPUS_ON_NODE")
    if slurm is not None and str(slurm).strip() != "":
        return max(1, int(slurm))
    n = os.cpu_count()
    return max(1, n if n is not None else 1)


def _strip_project_path(s):
    return (s or "").strip().rstrip("/\\")


load_dotenv()
docker_repo = "docker://mrcieu/genehackman"
user = os.getenv('USER')
input_file = "input.yaml"
start_time = datetime.now()

_project = _strip_project_path(os.getenv("PROJECT_DIR", ""))

DATA_DIR = format_dir_string(os.path.join(_project, "data"))
RESULTS_DIR = format_dir_string(os.path.join(_project, "results"))

_pipeline_log_override = os.getenv("PIPELINE_LOG_DIR")
_local_flag = os.getenv("GENEHACKMAN_LOCAL", "").strip().lower() in ("1", "true", "yes")
if _pipeline_log_override:
    pipeline_log_directory = format_dir_string(_pipeline_log_override)
elif _local_flag or not (DATA_DIR and os.path.isdir(DATA_DIR)):
    pipeline_log_directory = format_dir_string(os.path.join(os.getcwd(), "snakemake_logs"))
else:
    pipeline_log_directory = format_dir_string(f"{DATA_DIR}/snakemake_logs")

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

def _genehackman_package_version():
    description = Path(__file__).resolve().parent.parent.parent / "DESCRIPTION"
    if not description.is_file():
        return None
    for line in description.read_text(encoding="utf-8").splitlines():
        if line.startswith("Version:"):
            version = line.split(":", 1)[1].strip()
            return version or None
    return None


_docker_version_env = (os.getenv("DOCKER_VERSION") or "").strip()
DOCKER_VERSION = _docker_version_env or _genehackman_package_version()

# Single Snakemake output path for successful finemap runs (written last by finemap_gwas).
FINEMAP_COMPLETE_TXT_PATTERN = RESULTS_DIR + "finemap/{prefix}/finemap_complete.txt"

PIPELINE_DATA_DIR = format_dir_string(os.getenv('PIPELINE_DATA_DIR'))
QTL_DATA_DIR = os.getenv("QTL_DATA_DIR", "").strip()

LDSC_DIR = format_dir_string(PIPELINE_DATA_DIR + "/LDSCORE/b37_dbsnp156")
GENOMIC_DATA_DIR = format_dir_string(PIPELINE_DATA_DIR + "/genomic_data")
THOUSAND_GENOMES_DIR = format_dir_string(GENOMIC_DATA_DIR + "1000genomes/b37_dbsnp156")
QTL_DIRECTORY = format_dir_string(QTL_DATA_DIR)

AVAILABLE_CPUS = _available_cpus()
