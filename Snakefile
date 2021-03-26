# snakemake --jobs 5000 --latency-wait 90 --cluster-config cluster.json --cluster
# 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.14.0")


##### load config and sample sheets #####
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"],
    drop=False)
# enforce str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
validate(units, schema="schemas/units.schema.yaml")


##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.csv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg",
        "qc/multiqc_report.html"

##### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
