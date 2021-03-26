log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggplot2")
library("cowplot")
library("limma")
library("ggrepel")

dds <- readRDS(snakemake@input[[1]])
pca_color <- snakemake@params[['color']]
pca_fill <- snakemake@params[['fill']]


if (all(c(pca_color, pca_fill) != "")) {
    intgroup <- c(pca_color, pca_fill)
} else if (pca_color != "") {
    intgroup <- pca_color
} else if (pca_fill != "") {
    intgroup <- pca_fill
} else {
    stop("At least one of fill or color have to be specified")
}

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
#assay(counts) <- limma::removeBatchEffect(assay(counts), counts$individual)
pcaData <- plotPCA(counts, intgroup = intgroup, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


if (all(c(pca_color, pca_fill) != "")) {
    color_sym = sym(pca_color)
    fill_sym = sym(pca_fill)
    p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!pca_color,
                                                          fill=!!pca_fill))
        p <- p + geom_point(pch=22,  size=3, stroke=2)
        p <- p + geom_text(aes(label=name), check_overlap=TRUE)
} else {
        intgroup_sym = sym(intgroup)
    p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!intgroup_sym))
        p <- p + geom_point() +
                    scale_color_brewer(type="qual", palette="Dark2")
                    p <- p + geom_text(aes(label=name), check_overlap=TRUE)
}

p <- p +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.justification = 0)
ggsave(plot=p, height=4.5, width=7.5, filename = "results/pca.pdf")
ggsave(plot=p, height=4.5, width=7.5, filename = "results/pca.svg")


