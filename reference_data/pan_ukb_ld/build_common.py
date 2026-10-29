"""Shared helpers for the Pan-UKB LD pack build scripts (paths, logging, source access, Hail)."""

import concurrent.futures
import json
import logging
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
sys.path.insert(0, HERE)

AUTOSOMES = [str(c) for c in range(1, 23)]


def env(name, default=None):
    value = os.environ.get(name)
    return default if value is None or value == "" else value


def ld_settings():
    """The settings the LD of a chromosome depends on, in the same form as run_build.sh's settings_for ld.

    Written into state/<POP>.ld.chr<C>.done, so a chromosome built with other settings (e.g. phase 0 at
    another radius) is recognised as stale instead of being reused."""
    return (f"maf={env('MAF_MIN', '0.01')} radius={env('RADIUS_KB', '1000')} "
            f"deadzone={env('DEADZONE', '0.02')} chunk={env('CHUNK_ROWS', '2048')}")


class Paths:
    """Directory layout under WORK_DIR (see README.md)."""

    def __init__(self, work_dir=None):
        self.work = os.path.abspath(work_dir or env("WORK_DIR", "/local-scratch/data/pan_ukb_ld"))
        self.source = os.path.join(self.work, "source")
        self.raw = os.path.join(self.work, "raw")
        self.chunks = os.path.join(self.work, "work")
        self.pack = os.path.join(self.work, "pack")
        self.phase0 = os.path.join(self.work, "phase0")
        self.state = os.path.join(self.work, "state")
        self.logs = os.path.join(self.work, "logs")
        self.tmp = os.path.join(self.work, "tmp")

    def make(self):
        for path in (self.work, self.source, self.raw, self.chunks, self.pack, self.phase0,
                     self.state, self.logs, self.tmp):
            os.makedirs(path, exist_ok=True)
        return self

    def source_pop(self, pop):
        return os.path.join(self.source, pop)

    def raw_pop(self, pop):
        path = os.path.join(self.raw, pop)
        os.makedirs(path, exist_ok=True)
        return path

    def variant_ht(self, pop):
        return os.path.join(self.source_pop(pop), f"UKBB.{pop}.ldadj.variant.ht")

    def ldscore_ht(self, pop):
        return os.path.join(self.source_pop(pop), f"UKBB.{pop}.ldscore.ht")

    def block_matrix(self, pop):
        return os.path.join(self.source_pop(pop), f"UKBB.{pop}.ldadj.bm")


def source_names(pop):
    """Object prefixes (relative to SOURCE_URI) of the Pan-UKB files for one population."""
    return {"variant_ht": f"UKBB.{pop}.ldadj.variant.ht",
            "ldscore_ht": f"UKBB.{pop}.ldscore.ht",
            "block_matrix": f"UKBB.{pop}.ldadj.bm"}


def get_logger(name, log_file=None):
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handler = logging.FileHandler(log_file)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def read_json(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)


def free_bytes(path):
    os.makedirs(path, exist_ok=True)
    return shutil.disk_usage(path).free


# --------------------------------------------------------------------------- source stores

class LocalStore:
    """A directory laid out like the bucket (used for tests and for pre-downloaded mirrors)."""

    def __init__(self, root):
        self.root = root

    def list(self, prefix):
        base = os.path.join(self.root, prefix)
        out = []
        if os.path.isfile(base):
            return [(prefix, os.path.getsize(base))]
        for directory, _, files in os.walk(base):
            for name in files:
                full = os.path.join(directory, name)
                out.append((os.path.relpath(full, self.root).replace(os.sep, "/"), os.path.getsize(full)))
        return sorted(out)

    def download(self, key, dest):
        shutil.copyfile(os.path.join(self.root, key), dest)

    def read_bytes(self, key):
        with open(os.path.join(self.root, key), "rb") as handle:
            return handle.read()


class S3Store:
    """Anonymous (no-sign-request) access to a public S3 bucket via boto3."""

    def __init__(self, bucket, prefix):
        import boto3
        from botocore import UNSIGNED
        from botocore.config import Config
        self.bucket = bucket
        self.prefix = prefix.strip("/")
        self.client = boto3.client("s3", config=Config(signature_version=UNSIGNED, max_pool_connections=64,
                                                       retries={"max_attempts": 20, "mode": "adaptive"}))

    def _key(self, key):
        return f"{self.prefix}/{key}" if self.prefix else key

    def list(self, prefix):
        full = self._key(prefix)
        paginator = self.client.get_paginator("list_objects_v2")
        out = []
        for page in paginator.paginate(Bucket=self.bucket, Prefix=full):
            for item in page.get("Contents", []):
                rel = item["Key"][len(self.prefix) + 1:] if self.prefix else item["Key"]
                out.append((rel, int(item["Size"])))
        return sorted(out)

    def download(self, key, dest):
        from boto3.s3.transfer import TransferConfig
        self.client.download_file(self.bucket, self._key(key), dest,
                                  Config=TransferConfig(max_concurrency=4, multipart_chunksize=64 * 1024 ** 2))

    def read_bytes(self, key):
        return self.client.get_object(Bucket=self.bucket, Key=self._key(key))["Body"].read()


def open_store(uri=None):
    uri = uri or env("SOURCE_URI", "s3://pan-ukb-us-east-1/ld_release")
    if uri.startswith("s3://"):
        bucket, _, prefix = uri[len("s3://"):].partition("/")
        return S3Store(bucket, prefix)
    if uri.startswith("file://"):
        uri = uri[len("file://"):]
    return LocalStore(uri)


def download_files(store, items, dest_root, strip_prefix, logger, threads=16, retries=5):
    """Download (key, size) items under dest_root, skipping files already complete; resumable."""
    todo = []
    for key, size in items:
        rel = key[len(strip_prefix):].lstrip("/") if strip_prefix else key
        dest = os.path.join(dest_root, rel)
        if os.path.exists(dest) and os.path.getsize(dest) == size:
            continue
        todo.append((key, size, dest))
    total = sum(size for _, size, _ in todo)
    logger.info(f"downloading {len(todo)} files ({total / 1e9:.1f} GB); {len(items) - len(todo)} already present")
    if not todo:
        return 0

    def fetch(item):
        key, size, dest = item
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        tmp = dest + ".part"
        for attempt in range(1, retries + 1):
            try:
                store.download(key, tmp)
                if os.path.getsize(tmp) != size:
                    raise IOError(f"size mismatch for {key}: {os.path.getsize(tmp)} != {size}")
                os.replace(tmp, dest)
                return size
            except Exception as error:  # noqa: BLE001 - retry any transfer error
                if attempt == retries:
                    raise
                logger.warning(f"retry {attempt} for {key}: {error}")
                time.sleep(min(60, 5 * attempt))
        return 0

    done_bytes, done_files, started = 0, 0, time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        for size in pool.map(fetch, todo):
            done_bytes += size
            done_files += 1
            if done_files % 200 == 0 or done_files == len(todo):
                rate = done_bytes / max(1e-9, time.time() - started) / 1e6
                logger.info(f"  {done_files}/{len(todo)} files, {done_bytes / 1e9:.1f}/{total / 1e9:.1f} GB "
                            f"({rate:.0f} MB/s)")
    return total


# --------------------------------------------------------------------------- Hail

def init_hail(paths, log_name):
    import hail as hl
    cores = env("HAIL_CORES", "8")
    memory = env("HAIL_DRIVER_MEMORY", "32g")
    os.environ.setdefault("PYSPARK_SUBMIT_ARGS", f"--driver-memory {memory} pyspark-shell")
    os.makedirs(paths.tmp, exist_ok=True)
    hl.init(master=f"local[{cores}]", quiet=True, tmp_dir=paths.tmp, local_tmpdir=paths.tmp,
            log=os.path.join(paths.logs, f"hail_{log_name}_{time.strftime('%Y%m%d_%H%M%S')}.log"),
            spark_conf={"spark.driver.memory": memory, "spark.local.dir": paths.tmp,
                        "spark.driver.maxResultSize": "0"},
            default_reference="GRCh37")
    return hl
