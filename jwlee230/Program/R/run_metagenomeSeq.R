rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Read count data"),
                   make_option(c("-m", "--metadata"), type="character", default=NULL, help="Metadata"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input, metadata, output)
{
    library(metagenomeSeq)
    library(biomformat)
    library(dplyr)

    header <- read.table(input, header=FALSE, sep="\t", comment.char="", as.is=TRUE, nrows=1)
    b_data <- read.table(input, header=FALSE, sep="\t", comment.char="", as.is=TRUE, skip=1)
    colnames(b_data) <- unlist(header)
    rownames(b_data) <- b_data[, 1]
    print(head(b_data))

    header <- read.table(metadata, header=FALSE, sep="\t", comment.char="", as.is=TRUE, nrows=1)
    clin <- read.table(metadata, header=FALSE, sep="\t", comment.char="", as.is=TRUE, skip=1)
    colnames(clin) <- unlist(header)
    rownames(clin) <- clin[, 1]
    print(head(clin))

    ord <- intersect(colnames(b_data), rownames(clin))
    clin <- clin[ord, ]
    b_data <- b_data[, ord]

    obj <- newMRexperiment(b_data, phenoData=AnnotatedDataFrame(clin))
    obj <- cumNorm(obj, p=cumNormStatFast(obj))
    mod <- model.matrix(~Comparing, data=pData(obj))
    # res <- MRcoefs(fitZig(obj, mod, control=zigControl(maxit=99), useMixedModel=TRUE), number=100000, numberEff=TRUE, group=3, file=output)
    res <- MRcoefs(fitFeatureModel(obj, mod), number=100000, numberEff=TRUE, group=3, file=output)
}

if (length(opt) == 4)
{
    main(opt$input, opt$metadata, opt$output)
} else {
    print_help(opt_parser)
}
