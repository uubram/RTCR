# Description: parses RTCR output (e.g. "results.tsv") for the R tcR package.
#
# tcR website: http://imminfo.github.io/tcr/

parse.rtcr <- function(.filename)
{
    f <- gzfile(.filename)
    records <- strsplit(tail(readLines(f), n = -1), split = "\t", fixed = T)
    close(f)
    df <- data.frame(Umi.count = NA, Umi.proportion = NA,
                     Read.count = as.integer(sapply(records, function(x)x[1])),
                     CDR3.nucleotide.sequence = sapply(records,
                                                       function(x)x[5]),
                     CDR3.amino.acid.sequence = sapply(records,
                                                       function(x)x[2]),
                     V.gene = sapply(records,
                                     function(x)strsplit(x[3],"*",T)[[1]][1]),
                     D.gene = "",
                     J.gene = sapply(records,
                                     function(x)strsplit(x[4],"*",T)[[1]][1]),
                     V.end = sapply(records, function(x)x[6]),
                     J.start = sapply(records, function(x)x[7]),
                     D5.end = -1,
                     D3.end = -1,
                     VD.insertions = -1,
                     DJ.insertions = -1,
                     Total.insertions = -1,
                     stringsAsFactors = F)
    df$Read.proportion <- df$Read.count / sum(df$Read.count)
    df <- df[, c("Umi.count", "Umi.proportion", "Read.count",
                 "Read.proportion", "CDR3.nucleotide.sequence",
                 "CDR3.amino.acid.sequence", "V.gene", "J.gene", "D.gene",
                 "V.end", "J.start", "D5.end", "D3.end", "VD.insertions",
                 "DJ.insertions", "Total.insertions")]
    cls <- c("as.integer", "as.numeric", "as.integer", "as.numeric",
             "as.character", "as.character", "as.character", "as.character",
             "as.character", "as.integer", "as.integer", "as.integer",
             "as.integer", "as.integer", "as.integer", "as.integer")
    for (i in 1:ncol(df)){
        df[[i]] <- do.call(cls[i], list(df[[i]]))
    }
    df
}
