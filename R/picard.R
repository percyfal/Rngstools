setdec <- function(lines) {
    res = "."
    if (length(grep(",", lines)) > 0) {
        res=","
    }
    res
}


read.dup_metrics <- function(infile, pct.mult=TRUE, pct.char=list("PERCENT"), ...) {
    lines <- readLines(infile)
    res <- list()
    dec <- setdec(lines[7:8])
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", na.strings=c("?"), dec=dec, ...))
    if (pct.mult) {
        i <- unique(do.call("c", lapply(pct.char, function(x) {which(grepl(x, colnames(res$metrics)))})))
        res$metrics[,i] <- res$metrics[,i] * 100
    }
    res
}

read.align_metrics <- function(infile, pct.mult=TRUE, pct.char=list("PCT", "STRAND_BALANCE", "PF_HQ_ERROR_RATE"), ...) {
    lines <- readLines(infile)
    res <- list()
    dec <- setdec(lines[7:8])
    res$metrics <- rbind(read.table(textConnection(lines[7:10]), header=TRUE, fill=TRUE, sep="\t", na.strings=c("?"), dec=dec, ...))
    res$metrics <- cbind(CATEGORY=res$metrics[,1], as.data.frame(lapply(res$metrics[,-1], function(x){if (mode(x) != "numeric") x <- as.numeric(x); x})))
    res$metrics$CATEGORY <- factor(res$metrics$CATEGORY, levels=c("FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"))
    if (pct.mult) {
        i <- unique(do.call("c", lapply(pct.char, function(x) {which(grepl(x, colnames(res$metrics)))})))
        res$metrics[,i] <- res$metrics[,i] * 100
    }
    res
}

read.insert_metrics <- function(infile, ...) {
    lines <- readLines(infile)
    res <- list()
    dec <- setdec(lines[7:8])
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", na.strings=c("?"), dec=dec, ...))
    res$histogram <- rbind(read.table(textConnection(lines[11:length(lines)]), header=TRUE, fill=TRUE, sep="\t",...))
    res
}

read.hs_metrics <- function(infile, pct.mult=TRUE, pct.char=list("PCT", "ON_BAIT_VS_SELECTED"), ...) {
    lines <- readLines(infile)
    res <- list()
    dec = setdec(lines[7:8])
    res$metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", na.strings=c("?"), dec=dec, ...))
    if (pct.mult) {
        i <- unique(do.call("c", lapply(pct.char, function(x) {which(grepl(x, colnames(res$metrics)))})))
        res$metrics[,i] <- res$metrics[,i] * 100
    }
    v159cols <- c("AT_DROPOUT","GC_DROPOUT","SAMPLE","LIBRARY","READ_GROUP")
    if (length(colnames(res$metrics))==35) {
        res$metrics[,v159cols] <- NA
    }
    res
}

