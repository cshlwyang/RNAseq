def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return {"fq1":units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna(),
                "fq2":units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()}
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return {"fq1": expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=1, **wildcards),
                    "fq2": expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=2, **wildcards)}
        # single end sample
        return {"fq1":"trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)}

def get_fileend(wildcards):
    if config["trimming"]["skip"]:
        fq1=units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    else:
        if not is_single_end(**wildcards):
            fq1=expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=1, **wildcards)
        else:
            fq1=expand("trimmed/{sample}-{unit}.fastq.gz", **wildcards),
    if fq1[0].endswith(".gz"):
        readcmd = "--readFilesCommand zcat"
    else:
        readcmd = ""
    return readcmd


rule align:
    input:
        unpack(get_fq)
    output:
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        readcmd=get_fileend,
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    conda:
        "../envs/align.yaml"
    threads: 4
    shell:
        """
        STAR \
        {params.extra} \
            --runThreadN {threads} \
            --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} {params.readcmd} \
            --outReadsUnmapped Fastq \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix star/{wildcards.sample}-{wildcards.unit}/
        """
