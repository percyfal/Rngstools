## seqcap.R

get_miseq_metrics <- function(dir, runinfo, pattern) {
    metrics <- list.files(dir, pattern=pattern)
    laneids <- as.integer(gsub("^([0-9]+)_.*$", "\\1", metrics))
    names(metrics) <- runinfo$description[match(laneids, runinfo$lane)]
    metrics
}

flowcellreport.miseq <- function(indir, outdir, run_info="run_info.yaml") {
    fc <- basename(indir)
    tmp <- yaml.load_file(file.path(indir, run_info))
    runinfo <- as.data.frame(do.call("rbind", tmp$details))
    dupmetrics <- get_miseq_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.dup_metrics")
    alnmetrics <- get_miseq_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.align_metrics")
    hsmetrics <- get_miseq_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.hs_metrics")
    insertmetrics <- get_miseq_metrics(file.path(indir, "Analysis"), runinfo, pattern="*.insert_metrics")

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


get_hiseq_metrics <- function(dir, runinfo, pattern, type=c("dup", "hs", "insert", "align")) {
    metrics <- list.files(dir, pattern=pattern)
    sampleinfo <- lapply(strsplit(metrics, "_"), function(x){paste(x[1], gsub("([0-9]+)-.*", "\\1", x[5]), sep="_")})
    names(metrics) <- sampleinfo
    type <- match.arg(type)
    metrics.res <- lapply(metrics, function(x)  { f=file.path(dir, x);
                                                   switch(type, dup = read.dup_metrics(f), hs=read.hs_metrics(f),
                                                          insert = read.insert_metrics(f), align=read.align_metrics(f))})
    ## samples <- rownames(runinfo)
    ## if (type == "align") {
    ##     samples <- paste(rep(samples, each=3), 1:3, sep=".")
    ## }
    ## missingsamples <- samples[!samples %in% names(res)]
    metrics.res
}


flowcellreport.hiseq <- function(analysisdir, flowcelldir, outdir, run_info="run_info.yaml") {
    fc <- basename(analysisdir)
    tmp <- yaml.load_file(file.path(flowcelldir, run_info))
    runinfo.tmp <- as.data.frame(do.call("rbind", lapply(tmp, function(x) {do.call("cbind", x)})))
    mpnames <- c("barcode_id", "barcode_type", "genomes_filter_out","name", "sample_prj", "sequence")
    runinfo <- cbind(runinfo.tmp, do.call("rbind", lapply(runinfo.tmp$multiplex, function(y) {y[mpnames]})))
    rownames(runinfo) <- paste(runinfo$lane, runinfo$barcode_id, sep="_")
    samples <- rownames(runinfo)
    dupmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern="*.dup_metrics", type="dup")
    alnmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern="*.align_metrics", type="align")
    hsmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern="*.hs_metrics", type="hs")
    insertmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern="*.insert_metrics", type="insert")

    dupmetrics.res.tab <- do.call("rbind", lapply(dupmetrics.res, function(x) {x$metrics}))
    if (length(setdiff(rownames(runinfo), rownames(dupmetrics.res.tab))) > 0) {
        dupmetrics.res.tab[rownames(runinfo[!rownames(runinfo) %in% rownames(dupmetrics.res.tab),]),] <- NA
    }
    dupmetrics.res.tab$SAMPLE <- rownames(dupmetrics.res.tab)
    dupmetrics.res.tab$lane <- do.call("c", runinfo[match(dupmetrics.res.tab$SAMPLE, rownames(runinfo)),]$lane)
    dupmetrics.res.tab$project <- do.call("c", runinfo[match(dupmetrics.res.tab$SAMPLE,rownames(runinfo)),]$sample_prj)

    pdf(file=file.path(outdir, "dup-summary.pdf"))
    print(stripplot(PERCENT_DUPLICATION ~ SAMPLE, data=dupmetrics.res.tab,  scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Sample", ylab="Percent duplication"))
    dev.off()
    pdf(file=file.path(outdir, "dup-summary-lane.pdf"))
    print(stripplot(PERCENT_DUPLICATION ~ lane, data=dupmetrics.res.tab,  groups=project, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Lane", ylab="Percent duplication"))
    dev.off()
    pdf(file=file.path(outdir, "dup-summary-project.pdf"))
    print(stripplot(PERCENT_DUPLICATION ~ project, data=dupmetrics.res.tab,  groups=lane, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Project", ylab="Percent duplication"))
    dev.off()

    ## alnmetrics
    alnmetrics.res.tab <- do.call("rbind", lapply(alnmetrics.res, function(x) {x$metrics}))
    if (length(setdiff(paste(samples, "1", sep="."), rownames(alnmetrics.res.tab))) > 0) {
        missingsamples <- paste(rep(rownames(runinfo)[!paste(rownames(runinfo), 1, sep=".") %in% rownames(alnmetrics.res.tab)], each =3), 1:3, sep=".")
        alnmetrics.res.tab[missingsamples,] <- NA
    }
    alnmetrics.res.tab <- cbind(alnmetrics.res.tab, SAMPLE=do.call("cbind", strsplit(rownames(alnmetrics.res.tab), "\\.") )[1,])
    alnmetrics.res.tab$lane <- do.call("c", runinfo[match(alnmetrics.res.tab$SAMPLE, rownames(runinfo)),]$lane)
    alnmetrics.res.tab$project <- do.call("c", runinfo[match(alnmetrics.res.tab$SAMPLE,rownames(runinfo)),]$sample_prj)

    pdf(file=file.path(outdir, "mapping-summary.pdf"))
    print(stripplot(TOTAL_READS/1e6 + PF_READS_ALIGNED/1e6  ~ CATEGORY | SAMPLE, data=alnmetrics.res.tab, auto.key=list(text=c("Reads", "Aligned")), scales=list(x=list(rot=45)), ylab="Reads (millions)", xlab="Category.", par.settings=simpleTheme(col=c("black","red"), pch=21)))
dev.off()

    pdf(file=file.path(outdir, "mapping-summary-by-category.pdf"))
    print(stripplot(PCT_PF_READS_ALIGNED ~ SAMPLE | CATEGORY,groups=SAMPLE, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    pdf(file=file.path(outdir, "mapping-summary-by-lane.pdf"))
    print(stripplot(PCT_PF_READS_ALIGNED ~ lane | CATEGORY, groups=project, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45), y=list(rot=45, at=NULL))), par.settings=simpleTheme(pch=19)))
    dev.off()

    pdf(file=file.path(outdir, "mapping-summary-by-project.pdf"))
    print(stripplot(PCT_PF_READS_ALIGNED ~ project | CATEGORY, groups=lane, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    insertmetrics.res.tab <- do.call("rbind", lapply(insertmetrics.res, function(x) {x$metrics}))
    if (length(setdiff(rownames(runinfo), rownames(insertmetrics.res.tab))) > 0) {
        insertmetrics.res.tab[rownames(runinfo[!rownames(runinfo) %in% rownames(insertmetrics.res.tab),]),] <- NA
    }
    insertmetrics.res.hist <- do.call("rbind", lapply(names(insertmetrics.res), function(x) {cbind(insertmetrics.res[[x]]$histogram[,c("insert_size", "fr_count")], sample=x)}))

    pdf(file=file.path(outdir, "insert-summary.pdf"))
    print(xyplot(fr_count ~ insert_size | sample, data=insertmetrics.res.hist, xlab="Insert size", ylab="Count", xlim=c(0,1000), main="Insert size distributions", type="l", lwd=2))
    dev.off()

    hsmetrics.res.tab <- do.call("rbind", lapply(hsmetrics.res, function(x) {x$metrics}))
    if (length(setdiff(rownames(runinfo), rownames(hsmetrics.res.tab))) > 0) {
        hsmetrics.res.tab[rownames(runinfo[!rownames(runinfo) %in% rownames(hsmetrics.res.tab),]),] <- NA
    }
    hsmetrics.res.tab <- cbind(SAMPLE=rownames(hsmetrics.res.tab), hsmetrics.res.tab)
    hsmetrics.res.tab$lane <- do.call("c", runinfo[match(hsmetrics.res.tab$SAMPLE, rownames(runinfo)),]$lane)
    hsmetrics.res.tab$project <- do.call("c", runinfo[match(hsmetrics.res.tab$SAMPLE,rownames(runinfo)),]$sample_prj)
    hsmetrics.res.tab <- cbind(hsmetrics.res.tab, PERCENT_ON_TARGET=hsmetrics.res.tab$FOLD_ENRICHMENT/(hsmetrics.res.tab$GENOME_SIZE / hsmetrics.res.tab$TARGET_TERRITORY)* 100)


    pdf(file=file.path(outdir, "hs-summary.pdf"))
    print(xyplot(FOLD_ENRICHMENT / (GENOME_SIZE/TARGET_TERRITORY) ~ SAMPLE, data=hsmetrics.res.tab, xlab="Sample", ylab="Percent on target", main="Percent on target", scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    pdf(file=file.path(outdir, "hs-summary-by-project.pdf"))
    print(stripplot(100*FOLD_ENRICHMENT / (GENOME_SIZE/TARGET_TERRITORY) ~ project, data=hsmetrics.res.tab, xlab="Project", ylab="Percent on target", main="Percent on target", scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    i <- c("ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X")
    hsmetrics.stack <- stack(hsmetrics.res.tab[,i])
    hsmetrics.stack <- cbind(hsmetrics.stack, SAMPLE=hsmetrics.res.tab$SAMPLE, project=hsmetrics.res.tab$project)
    hsmetrics.stack$ind <- factor(hsmetrics.stack$ind, levels=levels(hsmetrics.stack$ind)[c(5,3,1,2,4)])
    levels(hsmetrics.stack$ind) <- c("0X", "2X", "10X", "20X", "30X")
    pdf(file=file.path(outdir, "hs-coverage-by-sample.pdf"))
    print(stripplot(values ~ ind | SAMPLE, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
    dev.off()
    pdf(file=file.path(outdir, "hs-coverage-by-project.pdf"))
    print(stripplot(values ~ ind | project, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
    dev.off()
    pdf(file=file.path(outdir, "hs-coverage-by-project-bw.pdf"))
    print(bwplot(values ~ ind | project, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
    dev.off()

    res.df <- cbind(runinfo[samples,], dupmetrics.res.tab[samples,], insertmetrics.res.tab[samples,], hsmetrics.res.tab[samples,],
                    alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="PAIR",][paste(samples, "3", sep="."),],
                    alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="FIRST_OF_PAIR",][paste(samples, "1", sep="."),],
                    alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",][paste(samples, "2", sep="."),])
    res.df$lane <- as.integer(res.df$lane)
    res.df$barcode_id <- as.integer(res.df$barcode_id)
    res.df$name <- as.character(res.df$name)

    res.list <- list(dupmetrics=dupmetrics.res.tab,
                     insertmetrics=insertmetrics.res.tab,
                     hsmetrics=hsmetrics.res.tab,
                     alnmetrics=alnmetrics.res.tab,
                     alnmetrics_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="PAIR",],
                     alnmetrics_first_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="FIRST_OF_PAIR",],
                     alnmetrics_second_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",])


    ## Write a table with the most important information
    ## dupmetrics
    reportfile <- file.path(outdir, paste(fc, "-metrics.txt", sep=""))
    ## res.df
    reportcols <- c("SAMPLE", "project", "lane", "barcode_id", "name", "TOTAL_READS", "PERCENT_DUPLICATION", "MEAN_INSERT_SIZE", "GENOME_SIZE", "FOLD_ENRICHMENT", "PERCENT_ON_TARGET", "PCT_USABLE_BASES_ON_TARGET", "PCT_TARGET_BASES_10X", "PCT_PF_READS_ALIGNED")
    write.table(file=reportfile, res.df[,reportcols], sep="\t", row.names=FALSE)

    res.list

}


print.flowcellmetrics <- function(res, outcols=c(10, 8, 14, 34, 52, 53, 70)) {
    tmp <- res[,outcols]
    tmp$PCT_ON_TARGET <- res$FOLD_ENRICHMENT / (res$GENOME_SIZE/res$TARGET_TERRITORY) * 100
    tmp
}
