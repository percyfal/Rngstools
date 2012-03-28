## seqcap.R

get_metrics <- function(dir, runinfo, pattern) {
    metrics <- list.files(dir, pattern=pattern)
    laneids <- as.integer(gsub("^([0-9]+)_.*$", "\\1", metrics))
    names(metrics) <- runinfo$description[match(laneids, runinfo$lane)]
    metrics
}

flowcellreport.miseq <- function(indir, outdir, run_info="run_info.yaml") {
    fc <- basename(indir)
    tmp <- yaml.load_file(file.path(indir, run_info))
    runinfo <- as.data.frame(do.call("rbind", tmp$details))
    dupmetrics <- get_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.dup_metrics")
    alnmetrics <- get_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.align_metrics")
    hsmetrics <- get_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.hs_metrics")
    insertmetrics <- get_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.insert_metrics")

    dupmetrics.res <- lapply(dupmetrics, function(f) {read.dup_metrics(file.path(indir, "Analysis",f))})
    alnmetrics.res <- lapply(alnmetrics, function(f) {read.align_metrics(file.path(indir, "Analysis",f))})
    insertmetrics.res <- lapply(insertmetrics, function(f) {read.insert_metrics(file.path(indir, "Analysis",f))})
    hsmetrics.res <- lapply(hsmetrics, function(f) {read.hs_metrics(file.path(indir, "Analysis",f))})
    
    dupmetrics.res.tab <- do.call("rbind", lapply(dupmetrics.res, function(x) {x$metrics}))
    dupmetrics.res.tab$SAMPLE <- rownames(dupmetrics.res.tab)

    pdf(file=file.path(outdir, "dup-summary.pdf"))
    print(stripplot(PERCENT_DUPLICATION ~ SAMPLE, data=dupmetrics.res.tab,  scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Sample", ylab="Percent duplication"))
    dev.off()
    
    alnmetrics.res.tab <- do.call("rbind", lapply(alnmetrics.res, function(x) {x$metrics}))
    alnmetrics.res.tab <- cbind(alnmetrics.res.tab, NAME= rep(names(alnmetrics.res), each=3))

    pdf(file=file.path(outdir, "mapping-summary.pdf"))
    print(stripplot(TOTAL_READS/1e6 + PF_READS_ALIGNED/1e6  ~ CATEGORY | NAME, data=alnmetrics.res.tab, auto.key=list(text=c("Reads", "Aligned")), scales=list(x=list(rot=45)), ylab="Reads (millions)", xlab="Category.", par.settings=simpleTheme(col=c("black","red"), pch=21)))
dev.off()

    pdf(file=file.path(outdir, "mapping-summary-by-category.pdf"))
    print(stripplot(PCT_PF_READS_ALIGNED ~ NAME | CATEGORY,groups=NAME, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()
    
    insertmetrics.res.tab <- do.call("rbind", lapply(insertmetrics.res, function(x) {x$metrics}))
    insertmetrics.res.hist <- do.call("rbind", lapply(names(insertmetrics.res), function(x) {cbind(insertmetrics.res[[x]]$histogram, sample=x)}))

    pdf(file=file.path(outdir, "insert-summary.pdf"))
    print(xyplot(fr_count ~ insert_size | sample, data=insertmetrics.res.hist, xlab="Insert size", ylab="Count", main="Insert size distributions", type="l", lwd=2))
    dev.off()
    
    hsmetrics.res.tab <- do.call("rbind", lapply(hsmetrics.res, function(x) {x$metrics}))
    hsmetrics.res.tab <- cbind(SAMPLE=rownames(hsmetrics.res.tab), hsmetrics.res.tab)

    pdf(file=file.path(outdir, "hs-summary.pdf"))
    print(xyplot(FOLD_ENRICHMENT / (GENOME_SIZE/TARGET_TERRITORY) ~ SAMPLE, data=hsmetrics.res.tab, xlab="Sample", ylab="Percent on target", main="Percent on target", scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    res.df <- cbind(dupmetrics.res.tab, insertmetrics.res.tab, hsmetrics.res.tab,
                    as.matrix(alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="PAIR",], rownames=FALSE),
                    as.matrix(alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="FIRST_OF_PAIR",], rownames=FALSE),
                    as.matrix(alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",], rownames=FALSE))

    res.df
}

print.flowcellmetrics <- function(res, outcols=c(10, 8, 14, 34, 52, 53, 70)) {
    tmp <- res[,outcols]
    tmp$PCT_ON_TARGET <- res$FOLD_ENRICHMENT / (res$GENOME_SIZE/res$TARGET_TERRITORY) * 100
    tmp
}
