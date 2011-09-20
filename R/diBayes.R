## TODO move to data file
.titv <- as.factor(c("ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "ti", "tv",
                    "tv", "ti", "tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv",
                    "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "ti", "tv",
                    "tv", "ti", "tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv"))
.iub <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K")
.titv.df <- matrix(.titv, ncol=4)
colnames(.titv.df) <- .iub[1:4]
rownames(.titv.df) <- .iub

## All diBayes related methods are put here
setClass("diBayes", contains="RangedData", representation=representation(titv="character", funccode="factor", gene="character"))
setClass("diBayes", contains="RangedData")

## helper function to import gff to diBayes
read.diBayes.gff <- function(file) {
    gff <- import.gff(file)
    ## dbSNP function codes: Locus region = 1; Coding = 2; Coding-synon
    ## = 3; Coding-nonsynon = 4; mRNA-UTR = 5; Intron = 6; Splice-site
    ## = 7; Contig-reference = 8; Coding-exception = 9; NearGene-3 =
    ## 13; NearGene-5 = 15; Coding-nonsynonymous nonsense = 41;
    ## Coding-nonsynonymous missense = 42; Coding-nonsynonymous
    ## frameshift = 44; Coding-nonsynonymous cds-indel = 45; UTR-3 =
    ## 53; UTR-5 = 55; Splice-3 = 73; Splice-3 = 75.
    fc <- do.call("c", lapply(strsplit(gsub("[A-Za-z0-9\\-\\.]+[:\\-]", "", gff$functionCode), ","), function(x) { paste(as.vector(unique(sort(x))), collapse=",")}))
    fc[fc==""] <- "NA"
    gff$funccode <- factor(fc)
    gff$titv <- factor(mapply(function(x,y) {.titv.df[x,y]}, gff$genotype, gff$reference), levels=c("ti", "tv"))
    names(gff$titv) <- NULL
    ## Dangerous: assuming gene names have no commas or :
    gff$gene <- gsub(":[0-9]+$", "", gsub("^[0-9]+:", "", gsub(",.*", "", gff$functionCode)))
    gff.db <- new("diBayes", gff)
    gff.db
}

## filterVariants
setGeneric("filterVariants", function(object, cov=c(0,Inf), score=c(0, 1), genenames=NULL, maf.freq=0, only.novel=FALSE, only.dbsnp=FALSE) {standardGeneric("filterVariants")})

setMethod("filterVariants",
          signature=signature(
          object="diBayes"),
          function(object, cov=c(0,Inf), score=c(0, 1), genenames=NULL, maf.freq=0, only.novel=FALSE, only.dbsnp=FALSE) {
              i.all <- rep(FALSE, dim(object)[1])
              i.cov <- object$coverage >= cov[1] & object$coverage <= cov[2]
              i.score <- object$score >= score[1] & object$score <= score[2]
              i.maffreq <- object$novelAlleleCounts/(object$refAlleleCounts + object$novelAlleleCounts) > maf.freq
              i.names <- !i.all
              if (!is.null(genenames)) {
                  i.names <- i.all
                  i.names[do.call("c", lapply(genenames, grep, object$functionCode))] <- TRUE
              }
              i <- i.cov & i.score & i.names & i.maffreq
              if (only.novel)
                  i <- i & is.na(object$rsID)
              if (only.dbsnp)
                  i <- i & !is.na(object$rsID)
    new("diBayes", object[i,])
          }
          )

## plotQC
setGeneric("plotQC", function(object, formula = ~., range=c(0, max(object$coverage)), stepsize=1, titv.ratio=TRUE, return.tab=FALSE, ...) {
    standardGeneric("plotQC")})

setMethod("plotQC",
          signature=signature(
          object="diBayes"),
          function(object, formula = ~., range=c(0,max(object$coverage)), stepsize=1, titv.ratio=TRUE, return.tab=FALSE, ...) {
              ylab = "Freq"
              db = object[object$coverage >= min(range) & object$coverage <= max(range),]
              steps <- seq(min(range), max(range), by=stepsize)
              object$intervals <- cut(object$coverage, steps, include.lowest=TRUE)
              term.labels <- labels(terms.formula(formula))
              term.labels <- gsub("\\|", "+", term.labels)
              if (titv.ratio) {
                  term.labels <- gsub("[+ ]*titv", "", term.labels)
                  tab.df <- as.data.frame(xtabs(as.formula(paste("~ ", term.labels, "+ titv", sep="")), data=as.data.frame(db)))
                  n.row <- dim(tab.df)[1]/2
                  i.col <- match(c("titv", "Freq"), colnames(tab.df))
                  c.names <- colnames(tab.df[-i.col])
                  tab.df <- as.data.frame(cbind(tab.df[1:n.row,-i.col], Freq = tab.df$Freq[tab.df$titv=="ti"] / tab.df$Freq[tab.df$titv=="tv"]))
                  colnames(tab.df)[1:length(c.names)] <- c.names
                  newf <- as.formula(paste("Freq ~ ", gsub("\\|[ ]*titv", "", labels(terms.formula(formula)))))
                  ylab = "titv ratio"
              } else {
                  tab.df <- as.data.frame(xtabs(as.formula(paste("~ ", term.labels, sep="")), data=as.data.frame(db)))
                  newf <- as.formula(paste("Freq ~ ", labels(terms.formula(formula))))
              }
              if (return.tab)
                  tab.df
              else
                  xyplot(newf, data=tab.df, ylab=ylab, ...)
          }
          )

setGeneric("Barchart", function(object, formula, ...) {standardGeneric("Barchart")})

setMethod("Barchart",
          signature=signature(
          object="diBayes"),
          function(object, formula, ...) {
              tab.df <- xtabs(formula, data=object)
              barchart(tab.df, ...)
          }
          )

setMethod("summary",
          signature=signature(
          object="diBayes"),
          function(object) {
              obj.df <- as.data.frame(object)
              res <- list()
              labs <- c("hom", "het", "all")
              indbsnp <- factor(!is.na(object$rsID), levels=c(FALSE, TRUE))
              object$het <- factor(object$het, levels=c(0,1))
              dbsnp <- table(indbsnp, object$het)
              res$dbsnp <- list()
              res$dbsnp$conc <- c(prop.table(table(indbsnp, object$het),2)[2,],
                                  prop.table(margin.table(dbsnp, 1))[2])
              names(res$dbsnp$conc) <- labs
              res$dbsnp$count <- c(dbsnp[2,],
                                   margin.table(dbsnp,1)[2])
              names(res$dbsnp$count) <- labs
              res$basic <- xtabs(~ het, data=obj.df)
              res$basic <- c(res$basic, all=margin.table(res$basic))
              names(res$basic) <- labs
              genotype <-  factor(object$genotype, levels=.iub)
              res$snpchanges <- table(genotype, reference=factor(object$reference, levels=c("A", "C", "G", "T")) )
              res$titv <- xtabs(~ titv + het, data=obj.df)
              res$titv <- cbind(res$titv, all=margin.table(res$titv, 1))
              colnames(res$titv) <- labs
              ## Chromosome stats
              res$chrstats <- xtabs(~ space + het, data=obj.df)
              res$chrstats <- cbind(res$chrstats, all=margin.table(res$chrstats, 1))
              colnames(res$chrstats) <- labs
              res
          }
          )

setGeneric("funcCodes", function(object) {standardGeneric("funcCodes")})

setMethod("funcCodes",
          signature=signature(
          object="diBayes"),
          function(object) {
              i <- do.call("c", lapply(object$funccode, function(x) {y <- strsplit(as.character(x), ","); length(y[[1]])}))
              fc <- list(fc=do.call("c", lapply(object$funccode, function(x) {y <- strsplit(as.character(x), ","); y[[1]]})),
                         j = rep(1:length(i), i)
                         )
              fc
          }
          )

setGeneric("plotFunctionalCode", function(object, ...) {standardGeneric("plotFunctionalCode")})

setMethod("plotFunctionalCode",
          signature=signature(
          object="diBayes"),
          function(object, ...) {
              funccodes <- list("Locus region", "Coding","Coding-synon","Coding-nonsynon",
                                "mRNA-UTR", "Intron", "Splice-site",
                                "Contig-reference", "Contig-exception",
                                "NearGene-3", "NearGene-5",
                                "Coding-nonsynonymous nonsense", "Coding-nonsynonymous missense",
                                "Coding-nonsynonymous frameshift",
                                "Coding-nonsynonymous cds-indel", "UTR-3", "UTR-5", "Splice-3",
                                "Splice-5")
              names(funccodes) <- as.character(c(1,2,3,4,5,6,7,8,9,13,15,41,42,44,45,53,55,73,75))
              fc <- funcCodes(object)
              fc.tab <- table(fc$fc)
              i <- match(names(fc.tab), names(funccodes))
              labs <- do.call("rbind", funccodes[i[!is.na(i)]])
              labs <- paste(rownames(labs), labs)
              barplot(fc.tab, horiz=TRUE, legend.text=labs, ...)
          }
          )

##setGeneric("print", function(object, which="all") {standardGeneric("print")})
setMethod("print",
          signature=signature(
          object="diBayes"),
          function(object, which="all") {
              objsum <- summary(object)
              res <- rbind(all=objsum$basic, objsum$titv, dbsnp=objsum$dbsnp$count)
              res
          }
          )

