## TODO move to data file
.titv <- as.factor(c("ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "ti", "tv",
                    "tv", "ti", "tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv",
                    "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "ti", "tv",
                    "tv", "ti", "tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv"))
.iub <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K")

## TODO: Make S4 generic function on a RangedData object
plotQCdiBayes <-
function(db, steps=c(1,5,100), breaks=c(50, 100, 10000), ...) {
    if (length(steps) != length(breaks))
        stop("length(steps) != length(breaks)")
    cr <- range(db$coverage)
    s <- c(hist(db$coverage[db$coverage >= min(cr) & db$coverage < breaks[1]], plot=FALSE, breaks=(breaks[1]-min(cr))/steps[1])$breaks,
           hist(db$coverage[db$coverage >= breaks[1] & db$coverage < breaks[2]], plot=FALSE, breaks=(breaks[2]-breaks[1])/steps[2])$breaks,
           hist(db$coverage[db$coverage >= breaks[2] & db$coverage <= max(c(cr, breaks))], plot=FALSE, breaks=(max(c(cr, breaks)) - breaks[2])/steps[3])$breaks)
    slab <- paste(s-1, s, sep="-")
    dbqc <- lapply(s, function(x) {db.filt <- diBayes_filter(db, cov=c(x-1, x)); if (dim(db.filt)[1]>0) summarize_diBayes(db.filt) else NA})
    names(dbqc) <- slab
    titv <- do.call("rbind", lapply(names(dbqc), function(x) {dbqc[[x]]; if (!is.na(dbqc[[x]])) cbind(stack(as.data.frame(dbqc[[x]]$titv)), x)}))
}

## Enable formula expression?
plot.diBayes <- function(res, type="chr", genenames=NULL, cond=NULL, ...) {
    type <- match.arg(type, c("chr", "gene"))
    if (type == "chr") {
        chrdata <- space(res)
        i <- !grepl("[a-zA-Z]", levels(chrdata))
        newl <- c(sort(as.integer(levels(chrdata)[i])), levels(chrdata)[!i])
        chrdata <- factor(chrdata, levels=newl)
        barchart(chrdata, xlab="Freq", ...)
    }
    else if (type == "gene") {
        if (!is.null(genenames)) {
            z <- do.call("c", lapply(genenames, function(x) {tmp <- grep(x, res$functionCode); y <- if(length(tmp)>1) length(tmp) else 0; names(y) <- x; y}))
            barchart(z, xlab="Freq", ylab="Genes", ...)
        } else {
            stop("please provide gene names")
        }
    }
}

summary.diBayes <-
function(db, fc.return=FALSE) {
    res <- list()
    labs <- c("all", "hom", "het")
    indbsnp <- factor(!is.na(db$rsID), levels=c(FALSE, TRUE))
    db$het <- factor(db$het, levels=c(0,1))
    dbsnp <- table(indbsnp, db$het)
    res$dbsnp <- list()
    res$dbsnp$conc <- c(prop.table(margin.table(dbsnp, 1))[2],
                        prop.table(table(indbsnp, db$het),2)[2,])
    names(res$dbsnp$conc) <- labs
    res$dbsnp$count <- c(margin.table(dbsnp,1)[2],
                         dbsnp[2,])
    names(res$dbsnp$count) <- labs

    res$basic <- c(sum(table(db$het)), table(db$het))
    names(res$basic) <- labs

    genotype <-  factor(db$genotype, levels=.iub)
    res$snpchanges <- table(genotype, reference=factor(db$reference, levels=c("A", "C", "G", "T")) )
    het <- rep(c(rep("hom", 4), rep("het", 6)), 4)
    res$titv <- tapply(res$snpchanges, list(.titv, het), sum)
    res$titv <- cbind(res$titv, all=margin.table(res$titv, 1))
    ## Chromosome stats
    chrdata <- space(db)
    i <- !is.na(as.integer(levels(chrdata)))
    newl <- c(sort(as.integer(levels(chrdata))[i]), levels(chrdata)[!i])
    res$chrstats <- cbind(table(factor(chrdata, levels=newl)), table(factor(chrdata, levels=newl), db$het))
    colnames(res$chrstats) <- c("all", "hom", "het")

    if (fc.return) {
        ## Summary based on functional code
        i <- do.call("c", lapply(db$funccode, function(x) {y <- strsplit(as.character(x), ","); length(y[[1]])}))
        res$fc <- list(fc = do.call("c", lapply(db$funccode, function(x) {y <- strsplit(as.character(x), ","); y[[1]]})),
                       j = rep(1:length(i), i)
                       )
        class(res$fc) <- c("funccode", "list")
    } else {
        res$fc <- NULL
    }
    class(res) <- c("diBayes", "list")
    res
}

plot.funccode <-
function(fc, ...) {
    funccodes <- list("Locus region", "Coding","Coding-synon","Coding-nonsynon",
                      "mRNA-UTR", "Intron", "Splice-site",
                      "Contig-reference", "Contig-exception",
                      "NearGene-3", "NearGene-5",
                      "Coding-nonsynonymous nonsense", "Coding-nonsynonymous missense",
                      "Coding-nonsynonymous frameshift",
                      "Coding-nonsynonymous cds-indel", "UTR-3", "UTR-5", "Splice-3",
                      "Splice-5")
    names(funccodes) <- as.character(c(1,2,3,4,5,6,7,8,9,13,15,41,42,44,45,53,55,73,75))
    fctab <- table(fc$fc)
    i <- match(names(fctab), names(funccodes))
    labs <- do.call("rbind", funccodes[i[!is.na(i)]])
    labs <- paste(rownames(labs), labs)
    #if (length(labs) < length(names(fctab))) labs <- c(labs, NA)
    barplot(fctab, horiz=TRUE, legend.text=labs, ...)
}

filter.diBayes <-
function(db, cov=c(0,Inf), score=c(0, 1), genenames=NULL, ...) {
    i.all <- rep(FALSE, dim(db)[1])
    i.cov <- db$coverage >= cov[1] & db$coverage <= cov[2]
    i.score <- db$score >= score[1] & db$score <= score[2]
    i.names <- !i.all
    if (!is.null(genenames)) {
        i.names <- i.all
        i.names[do.call("c", lapply(genenames, grep, db$functionCode))] <- TRUE
    }
    i <- i.cov & i.score & i.names
    db[i,]
}

