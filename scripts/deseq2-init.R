log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("AnnotationHub")
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

message("Reading counts")
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene",
                  check.names=FALSE)
message("Reading sample file")
coldata <- read.table(snakemake@params[["samples"]], header=TRUE,
                      row.names="sample", check.names=FALSE, sep="\t")
message("Getting experimental design")
design <- as.formula(snakemake@params[["design"]])

# colData and countData must have the same sample order
if (nrow(coldata) != ncol(cts)) {
    stop("Number of samples in sample sheet and number of samples in counts",
         "matrix is not the same")
}
cts <- cts[,match(rownames(coldata),colnames(cts))]
if (any(c("control", "Control", "CONTROL") %in% levels(coldata$condition))) {
    if ("control" %in% levels(coldata$condition)) {
        coldata$condition <- relevel(coldata$condition, "control" )
    } else if ("Control" %in% levels(coldata$condition)) {
        coldata$condition <- relevel(coldata$condition, "Control" )
    } else {
        coldata$condition <- relevel(coldata$condition, "CONTROL" )
    }
}
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=design)

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 1,]

# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

# Remove build number on ENS gene id
rownames(dds) <- gsub("\\.\\d*", "", rownames(dds))

# Annotate by gene names
hub <- AnnotationHub()
# query(hub,  c("GTF", "Ensembl", "Mus musculus")) "AH7799"
# query(hub,  c("GTF", "Ensembl", "Homo sapiens")) "AH69461"
if (snakemake@params[["annotationhub"]] == "mouse") {
    hubid <- "AH7799"
} else if(snakemake@params[["annotationhub"]] == "human") {
    hubid <- "AH69461"
} else {
    stop("No annotation hub specified for organism:",
         snakemake@params[["annotationhub"]])
}
anno <- hub[[hubid]]
genemap <- tibble(gene_id=anno$gene_id,
                  symbol=anno$gene_name) %>%
    distinct()

featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    select(symbol)
mcols(dds) <- DataFrame(mcols(dds), featureData)

saveRDS(dds, file=snakemake@output[[1]])
