read.dup_metrics <- function(infile, ...) {
    lines <- readLines(infile)
    res <- list()
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", ...))
    res
}

read.align_metrics <- function(infile, ...) {
    lines <- readLines(infile)
    res <- list()
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", ...))
    res
}

read.insert_metrics <- function(infile, ...) {
    lines <- readLines(infile)
    res <- list()
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", ...))
    res$histogram <- rbind(read.table(textConnection(lines[11:length(lines)]), header=TRUE, fill=TRUE, sep="\t", ...))
    res
}

read.hs_metrics <- function(infile, ...) {
    lines <- readLines(infile)
    res <- list()
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", ...))
    res
}
