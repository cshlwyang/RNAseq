def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab",
        unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
       "../envs/pandas.yaml"
    log:
        "logs/counts/count_matrix.log"
    script:
        "../scripts/count-matrix.py"


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        "deseq2/all.rds"
    params:
        samples=config["samples"],
        design=config["diffexp"]["design"],
        annotationhub=config["diffexp"]["organism"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.rds"
    output:
        pca_svg=report("results/pca.svg", "../report/pca.rst"),
        pca_pdf="results/pca.pdf"
    params:
        fill=config["pca"]["fill"],
        color=config["pca"]["color"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.csv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
        ma_pdf="results/diffexp/{contrast}.ma-plot.pdf",
        up="results/diffexp/deg-sig-up_{contrast}.csv",
        down="results/diffexp/deg-sig-down_{contrast}.csv"
    params:
        contrast=get_contrast,
        annotationhub=config["diffexp"]["organism"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

