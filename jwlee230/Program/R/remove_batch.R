rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input TSV file", metavar="character"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output TSV files", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input, output)
{
    library(limma)
    input_data <- read.csv(input, stringsAsFactors=FALSE, sep="\t", skip=1, check.names=FALSE)

    OTU_ID = input_data[1]
    taxonomy_ID = input_data[length(input_data)]
    print(head(taxonomy_ID))

    input_data <- input_data[2:(length(input_data)-1)]
    first_batch <- vapply(strsplit(colnames(input_data), "-"), `[`, 1, FUN.VALUE=character(1))
    print(head(input_data))
    print(first_batch)

    output_data <- removeBatchEffect(input_data, batch=first_batch)
    output_data[output_data < 0] <- 0

    output_data <- cbind(OTU_ID, output_data, taxonomy_ID)

    cat("# Constructed from limma package\n", file=output)
    write.table(output_data, file=output, quote=FALSE, sep="\t", row.names=FALSE, append=TRUE)
}

if (length(opt) == 3)
{
    main(opt$input, opt$output);
} else {
    print_help(opt_parser)
}
