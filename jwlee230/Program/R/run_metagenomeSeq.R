
rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Read count data"),
                   make_option(c("-t", "--taxonomy"), type="character", default=NULL, help="Taxonomy data"),
                   make_option(c("-m", "--metadata"), type="character", default=NULL, help="Metadata"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input, taxonomy, metadata, output)
{
    library(metagenomeSeq)
    library(biomformat)
    library(dplyr)

    b <- biom2MRexperiment(read_biom(input))
    b_data <- biom_data(read_biom(input))
    print(rownames(b_data))

    taxa <- as.data.frame(read.table(taxonomy, header=TRUE, sep="\t", comment.char="", as.is=TRUE))
    rownames(taxa) <- taxa[, 1]
    print(head(taxa))

    header <- read.table(metadata, header=FALSE, sep="\t", comment.char="", as.is=TRUE, nrows=1)
    clin <- read.table(metadata, header=FALSE, sep="\t", comment.char="", as.is=TRUE, skip=2)
    colnames(clin) <- unlist(header)
    rownames(clin) <- clin[, 1]
    clin["SimplePremature"] <- clin["Simple Premature"]
    print(head(clin))

    ord <- intersect(colnames(b_data), rownames(clin))
    clin <- clin[ord, ]
    b_data <- b_data[, ord]

    ord <- intersect(rownames(b_data), rownames(taxa))
    taxa <- taxa[ord, ]
    b_data <- b_data[ord, ]

    obj <- newMRexperiment(b_data, phenoData=AnnotatedDataFrame(clin), featureData=AnnotatedDataFrame(taxa))
    obj <- cumNorm(obj, p=0.5)
    mod <- model.matrix(~1 + SimplePremature, data=pData(obj))
    res <- MRcoefs(fitFeatureModel(obj, mod))
    print(res)

    write.table(as.data.frame(res), file=output, sep="\t")
}

if (length(opt) == 5)
{
    main(opt$input, opt$taxonomy, opt$metadata, opt$output)
} else {
    print_help(opt_parser)
}
