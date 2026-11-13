"""Where each step of a pipeline gets its LD (and reference allele frequencies) from.

There is one LD source for the whole pipeline, the top-level ``ld_reference``, which every step that needs LD uses:

    harmonise   strand check and EAF filling (standardisation)
    clumping    clumping
    finemap     fine-mapping (SuSiE and MultiSuSiE)
    locus_zoom  r² and the index variant in regional plots (colocalization)

A step can use a different source (``<step>: {ld_reference: ...}``), and so can a GWAS
(``ld_reference`` on the GWAS entry, for example the genotypes it was computed from).  For a step and a GWAS the most
specific setting wins: the GWAS's, then the step's, then the top-level one.  There is no default: a step that runs
with nothing set stops the pipeline.

A source is ``pan_ukb`` (the Pan-UK Biobank LD), ``1000genomes`` (the 1000 Genomes panel of the GWAS's ancestry) or
the PLINK prefix (.bed/.bim/.fam) of the user's own genotypes, which may contain ``{ancestry}`` and environment
variables.  LDSC is not one of these steps: it uses LD scores (``ldsc.ld_scores``), not LD.

Plain Python (no Snakemake), so it can be tested on its own (tests/python/test_ld_sources.py).
"""
import os

LD_KEYWORDS = {
    "pan_ukb": "Pan-UK Biobank LD",
    "1000genomes": "1000 Genomes panel",
}
LD_STEPS = {
    "harmonise": "strand check and EAF filling",
    "clumping": "clumping",
    "finemap": "fine-mapping",
    "locus_zoom": "regional plots (LocusZoom r²)",
}
PLINK_EXTENSIONS = (".bed", ".bim", ".fam")

CHOICES = ("  ld_reference: pan_ukb        Pan-UK Biobank LD (~420,000 EUR; smaller for other ancestries)\n"
           "  ld_reference: 1000genomes    1000 Genomes panel of the GWAS's ancestry (~500 people per super-population)\n"
           "  ld_reference: /path/prefix   your own genotypes in PLINK format (.bed/.bim/.fam), in the output build;\n"
           "                               the path may contain {ancestry}, e.g. /data/genotypes_{ancestry}")


class LDSourceError(ValueError):
    """A problem with the LD settings in the input YAML; the message is meant for the user."""


def clean(value):
    """A setting as a stripped string ('' when not set)."""
    return "" if value is None else str(value).strip()


def is_genotypes(source):
    """True for a PLINK prefix (anything that is set and isn't one of the LD_KEYWORDS)."""
    return bool(source) and source not in LD_KEYWORDS


def genotype_prefix(source, ancestry=""):
    """Absolute PLINK prefix, with {ancestry}, ~ and environment variables filled in."""
    path = source.replace("{ancestry}", ancestry or "")
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def missing_plink_files(prefix):
    return [prefix + ext for ext in PLINK_EXTENSIONS if not os.path.isfile(prefix + ext)]


def resolve(step, top=None, step_value=None, gwas_value=None, ancestry="", where=""):
    """The LD source of one step for one GWAS: a keyword, an absolute PLINK prefix, or '' when nothing is set.

    gwas_value wins over step_value, which wins over top.  where: the GWAS, for messages."""
    raw = clean(gwas_value) or clean(step_value) or clean(top)
    if not is_genotypes(raw):
        return raw
    if "{ancestry}" in raw and not clean(ancestry):
        raise LDSourceError(f"ld_reference {raw} for {LD_STEPS.get(step, step)}{where} contains {{ancestry}}, but the "
                            "GWAS has no ancestry. Set its ancestry, or give the full path.")
    return genotype_prefix(raw, clean(ancestry))


def check_genotypes(source, step, where=""):
    """Stop when a PLINK prefix hasn't got its .bed/.bim/.fam files."""
    if not is_genotypes(source):
        return
    missing = missing_plink_files(source)
    if missing:
        raise LDSourceError(
            f"ld_reference for {LD_STEPS.get(step, step)}{where} is read as the PLINK prefix of your own genotypes "
            f"(it isn't {' or '.join(LD_KEYWORDS)}), but {', '.join(missing)} not found.\n"
            "Give the prefix without .bed/.bim/.fam, in a folder the container can see (your home folder, "
            "PIPELINE_DATA_DIR, DATA_DIR or RESULTS_DIR), or use one of:\n" + CHOICES)


def not_set_message(step, where=""):
    return (f"No LD reference is set for {LD_STEPS.get(step, step)}{where}. Choose where LD comes from: the top-level "
            "ld_reference is used by every step (strand check and EAF filling, clumping, fine-mapping and regional "
            f"plots), and '{step}: {{ld_reference: ...}}' or ld_reference on a GWAS sets it for one step or one GWAS.\n"
            + CHOICES)


def ld_bfile_message(gwas_file):
    return (f"ld_bfile (GWAS {gwas_file}) has been replaced: give the PLINK prefix as ld_reference on the GWAS instead. "
            "It is then used by every step for that GWAS (strand check and EAF filling, clumping, fine-mapping and "
            "regional plots), unless a step sets its own ld_reference.")


def label(source):
    """How a source is shown in the log and the reports."""
    if source in LD_KEYWORDS:
        return f"{LD_KEYWORDS[source]} ({source})"
    return f"your genotypes ({source})"


def describe(sources, ldsc_scores=None):
    """Lines saying where each step gets LD.

    sources: {step: {gwas label: source}} for the steps the pipeline runs, in LD_STEPS order.
    ldsc_scores: the LD score set used by LDSC (compare_gwases), shown as the exception."""
    lines = []
    for step in LD_STEPS:
        per_gwas = sources.get(step)
        if not per_gwas:
            continue
        distinct = sorted(set(per_gwas.values()))
        if len(distinct) == 1:
            lines.append(f"{LD_STEPS[step]}: {label(distinct[0])}")
        else:
            lines.append(f"{LD_STEPS[step]}: " + "; ".join(f"{g}: {label(s)}" for g, s in per_gwas.items()))
    if ldsc_scores:
        lines.append(f"LDSC: {ldsc_scores} LD scores (ldsc.ld_scores; LDSC uses LD scores, not ld_reference)")
    return lines
