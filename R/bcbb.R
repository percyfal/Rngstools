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


get_hiseq_metrics <- function(indir, runinfo, pattern, type=c("dup", "hs", "insert", "align"), sampletable=NULL, ...) {
    metrics <- list.files(indir, pattern=pattern)
    sampleinfo <- lapply(strsplit(metrics, "_"), function(x){paste(x[1], gsub("([0-9]+)-.*", "\\1", x[5]), sep="_")})
    sampleinfo <- gsub("_NA", "", sampleinfo)
    if (!is.null(sampletable)) {
        sampleinfo <- paste(sampletable[match(sampleinfo, sampletable[,1]),1], sampletable[match(sampleinfo, sampletable[,1]),2], sep="_")
    }
    names(metrics) <- sampleinfo
    type <- match.arg(type)
    metrics.res <- lapply(metrics, function(x)  { f=file.path(indir, x);
                                                   switch(type, dup = read.dup_metrics(f, ...), hs=read.hs_metrics(f, ...),
                                                          insert = read.insert_metrics(f, ...), align=read.align_metrics(f, ...))})
    ## samples <- rownames(runinfo)
    ## if (type == "align") {
    ##     samples <- paste(rep(samples, each=3), 1:3, sep=".")
    ## }
    ## missingsamples <- samples[!samples %in% names(res)]
    metrics.res
}



print.flowcellmetrics <- function(res, outcols=c(10, 8, 14, 34, 52, 53, 70)) {
    tmp <- res[,outcols]
    tmp$PCT_ON_TARGET <- res$FOLD_ENRICHMENT / (res$GENOME_SIZE/res$TARGET_TERRITORY) * 100
    tmp
}

flowcellreport.hiseq <- function(analysisdir, flowcelldir, outdir, run_info="run_info.yaml", dupmetrics=list(run=TRUE, pattern="*.dup_metrics", dec="."), insertmetrics=list(run=TRUE, pattern="*.insert_metrics"), alnmetrics=list(run=TRUE, pattern="*.align_metrics"), hsmetrics=list(run=TRUE, pattern="*.hs_metrics"), gcmetrics=list(run=TRUE, pattern="*.gc_metrics")) {

    fc <- basename(analysisdir)
    tmp <- yaml.load_file(file.path(flowcelldir, run_info))
    runinfo.tmp <- as.data.frame(do.call("rbind", lapply(tmp, function(x) {do.call("cbind", x)})))
    mpnames <- c("barcode_id", "barcode_type", "genomes_filter_out","name", "sample_prj", "sequence")
    runinfo <- cbind(runinfo.tmp, do.call("rbind", lapply(runinfo.tmp$multiplex, function(y) {y[mpnames]})))
    rownames(runinfo) <- paste(runinfo$lane, runinfo$barcode_id, sep="_")
    samples <- rownames(runinfo)
    res.list <- list(dupmetrics=NA, insertmetrics=NA, alnmetrics=NA,
                     alnmetrics_pair=NA, alnmetrics_first_of_pair=NA,
                     alnmetrics_second_of_pair=NA)

    res.df <- NULL
    ## duplication metrics
    if (dupmetrics$run) {
        cat("Adding dupmetrics...\n")
        dupmetrics.res.tab <- getDupmetrics(analysisdir, runinfo, outdir, dupmetrics$pattern)
        res.list$dupmetrics <- dupmetrics.res.tab
        res.df <- dupmetrics.res.tab
        reportfile <- file.path(outdir, "dupmetrics.txt")
        write.table(file=reportfile, dupmetrics.res.tab, sep="\t", row.names=TRUE)
    }

    ## alnmetrics
    if (alnmetrics$run) {
        cat("Adding alnmetrics...\n")
        alnmetrics.res.tab <- getAlnmetrics(analysisdir, runinfo, outdir, pattern=alnmetrics$pattern)
        res.list$alnmetrics <- alnmetrics.res.tab
        res.list$alnmetrics_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="PAIR",]
        colnames(res.list$alnmetrics_pair) <- paste("aln", colnames(res.list$alnmetrics_pair), "3", sep="_")
        res.list$alnmetrics_first_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="FIRST_OF_PAIR",]
        colnames(res.list$alnmetrics_first_of_pair) <- paste("aln", colnames(res.list$alnmetrics_first_of_pair), "1", sep="_")
        res.list$alnmetrics_second_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",]
        colnames(res.list$alnmetrics_second_of_pair) <- paste("aln", colnames(res.list$alnmetrics_second_of_pair), "2", sep="_")
        res.df <- cbind(res.df,
                        res.list$alnmetrics_pair,
                        res.list$alnmetrics_first_of_pair,
                        res.list$alnmetrics_second_of_pair)
                        #alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",][paste(samples, "2", sep="."),])
        reportfile <- file.path(outdir, "alnmetrics.txt")
        write.table(file=reportfile, alnmetrics.res.tab, sep="\t", row.names=TRUE)
    }
    ## insertmetrics
    if (insertmetrics$run) {
        cat("Adding insertmetrics...\n")
        insertmetrics.res.tab <- getInsertmetrics(analysisdir, runinfo, outdir, pattern=insertmetrics$pattern)
        res.list$insertmetrics=insertmetrics.res.tab
        res.df <- cbind(res.df, insertmetrics.res.tab)
        reportfile <- file.path(outdir, "insertmetrics.txt")
        write.table(file=reportfile, insertmetrics.res.tab, sep="\t", row.names=TRUE)
    }

    ## hsmetrics
    if (hsmetrics$run) {
        cat("Adding hsmetrics...\n")
        hsmetrics.res.tab <- getHsmetrics(analysisdir, runinfo,  outdir, pattern=hsmetrics$pattern)
        res.list$hsmetrics=hsmetrics.res.tab
        res.df <- cbind(res.df, hsmetrics.res.tab)
        reportfile <- file.path(outdir, "hsmetrics.txt")
        write.table(file=reportfile, hsmetrics.res.tab, sep="\t", row.names=TRUE)
    }
    ##res.df <- cbind(res.df, runinfo[match(rownames(res.df), rownames(runinfo)),])
    res.df$lane <- as.integer(res.df$lane)

    ## Write a table with the most important information
    ## dupmetrics
    reportfile <- file.path(outdir, paste(fc, "-metrics.txt", sep=""))
    reportcols <- c("SAMPLE", "project", "lane", "TOTAL_READS", "PERCENT_DUPLICATION", "MEAN_INSERT_SIZE", "GENOME_SIZE", "FOLD_ENRICHMENT", "PERCENT_ON_TARGET", "PCT_USABLE_BASES_ON_TARGET", "PCT_TARGET_BASES_10X", "aln_PCT_PF_READS_ALIGNED_1")
    write.table(file=reportfile, res.df[,reportcols], sep="\t", row.names=FALSE)
    list(res.list, res.df)
}

projectreport <- function(analysisdir, run_info, outdir, dupmetrics=list(run=TRUE, pattern="*.dup_metrics"), insertmetrics=list(run=TRUE, pattern="*.insert_metrics"), alnmetrics=list(run=TRUE, pattern="*.align_metrics"), hsmetrics=list(run=TRUE, pattern="*.hs_metrics"), gcmetrics=list(run=TRUE, pattern="*.gc_metrics")) {
    tmp <- yaml.load_file(run_info)
    if (names(tmp) == "details") {
        tmp <- tmp$details
    }
    runinfo <- as.data.frame(do.call("rbind", lapply(tmp, function(x) {y <- as.data.frame(rbind(x)); rownames(y) <- paste(x$lane, x$description, sep="_"); y})))
    fcinfo <- do.call("rbind", lapply(runinfo$files, function(x) {m <- regexec("(^[0-9]+).*_([0-9]{6})_([0-9A-Z]+)_.*", basename(x)[[1]]);y <- regmatches(basename(x)[[1]], m); y[[1]][2:4] }))
    colnames(fcinfo) <- c("fc_lane", "date", "flowcell")
    runinfo <- cbind(runinfo, fcinfo)
    samples <- runinfo$description
    res.list <- list(dupmetrics=NA, insertmetrics=NA, alnmetrics=NA,
                     alnmetrics_pair=NA, alnmetrics_first_of_pair=NA,
                     alnmetrics_second_of_pair=NA)
    res.df <- NULL
    ## duplication metrics
    if (dupmetrics$run) {
        cat("Adding dupmetrics...\n")
        dupmetrics.res.tab <- getDupmetrics(analysisdir, runinfo, outdir, dupmetrics$pattern)
        res.list$dupmetrics <- dupmetrics.res.tab
        res.df <- dupmetrics.res.tab
        reportfile <- file.path(outdir, "dupmetrics.txt")
        write.table(file=reportfile, dupmetrics.res.tab, sep="\t", row.names=FALSE)
    }
    ## alnmetrics
    if (alnmetrics$run) {
        cat("Adding alnmetrics...\n")
        alnmetrics.res.tab <- getAlnmetrics(analysisdir, runinfo, outdir, pattern=alnmetrics$pattern)
        res.list$alnmetrics <- alnmetrics.res.tab
        res.list$alnmetrics_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="PAIR",]
        colnames(res.list$alnmetrics_pair) <- paste("aln", colnames(res.list$alnmetrics_pair), "3", sep="_")
        res.list$alnmetrics_first_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="FIRST_OF_PAIR",]
        colnames(res.list$alnmetrics_first_of_pair) <- paste("aln", colnames(res.list$alnmetrics_first_of_pair), "1", sep="_")
        res.list$alnmetrics_second_of_pair=alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",]
        colnames(res.list$alnmetrics_second_of_pair) <- paste("aln", colnames(res.list$alnmetrics_second_of_pair), "2", sep="_")
        res.df <- cbind(res.df,
                        res.list$alnmetrics_pair,
                        res.list$alnmetrics_first_of_pair,
                        res.list$alnmetrics_second_of_pair)
                        #alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY=="SECOND_OF_PAIR",][paste(samples, "2", sep="."),])
        reportfile <- file.path(outdir, "alnmetrics.txt")
        write.table(file=reportfile, alnmetrics.res.tab, sep="\t", row.names=FALSE)
    }

    ## insertmetrics
    if (insertmetrics$run) {
        cat("Adding insertmetrics...\n")
        insertmetrics.res.tab <- getInsertmetrics(analysisdir, runinfo, outdir, pattern=insertmetrics$pattern)
        res.list$insertmetrics=insertmetrics.res.tab
        res.df <- cbind(res.df, insertmetrics.res.tab)
        reportfile <- file.path(outdir, "insertmetrics.txt")
        write.table(file=reportfile, insertmetrics.res.tab, sep="\t", row.names=FALSE)
    }

    ## hsmetrics
    if (hsmetrics$run) {
        cat("Adding hsmetrics...\n")
        hsmetrics.res.tab <- getHsmetrics(analysisdir, runinfo,  outdir, pattern=hsmetrics$pattern)
        res.list$hsmetrics=hsmetrics.res.tab
        res.df <- cbind(res.df, hsmetrics.res.tab)
        reportfile <- file.path(outdir, "hsmetrics.txt")
        write.table(file=reportfile, hsmetrics.res.tab, sep="\t", row.names=FALSE)
    }

    res.df <- cbind(res.df, runinfo[match(rownames(res.df), rownames(runinfo)),])
    res.df$lane <- as.integer(res.df$lane)

    list(res.list, res.df)
}



getDupmetrics <- function(analysisdir, runinfo, outdir, pattern="*.dup_metrics") {
    sample_prj <- FALSE
    if ("sample_prj" %in% names(runinfo)) {
        dupmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="dup")
        sample_prj <- TRUE
    } else {
        dupmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="dup", sampletable=runinfo[,c("lane", "description")])
    }
    dupmetrics.res.tab <- do.call("rbind", lapply(dupmetrics.res, function(x) {x$metrics}))
    sampleinfo <- getLaneAndDescription(dupmetrics.res.tab, concatenate=sample_prj)

    pdf(file=file.path(outdir, "dup-summary.pdf"))
    print(stripplot(PERCENT_DUPLICATION ~ sampleinfo$description, data=dupmetrics.res.tab,  scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Sample", ylab="Percent duplication"))
    dev.off()

    if ("sample_prj" %in% names(runinfo)) {
        dupmetrics.res.tab$SAMPLE <- rownames(dupmetrics.res.tab)
        dupmetrics.res.tab$lane <- do.call("c", runinfo[match(dupmetrics.res.tab$SAMPLE, rownames(runinfo)),]$lane)
        dupmetrics.res.tab$project <- do.call("c", runinfo[match(dupmetrics.res.tab$SAMPLE,rownames(runinfo)),]$sample_prj)
        pdf(file=file.path(outdir, "dup-summary-lane.pdf"))
        print(stripplot(PERCENT_DUPLICATION ~ lane, data=dupmetrics.res.tab,  groups=project, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Lane", ylab="Percent duplication"))
        dev.off()
        pdf(file=file.path(outdir, "dup-summary-project.pdf"))
        print(stripplot(PERCENT_DUPLICATION ~ project, data=dupmetrics.res.tab,  groups=lane, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),  ylim=c(0,100), xlab="Project", ylab="Percent duplication"))
        dev.off()
    }

    dupmetrics.res.tab
}

getAlnmetrics <- function(analysisdir, runinfo, outdir, pattern="*.align_metrics") {
    sample_prj <- FALSE
    if ("sample_prj" %in% names(runinfo)) {
        alnmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="align")
        sample_prj <- TRUE
    } else {
        alnmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="align", sampletable=runinfo[,c("lane", "description")])
    }
    alnmetrics.res.tab <- do.call("rbind", lapply(alnmetrics.res, function(x) {x$metrics}))
    sampleinfo <- getLaneAndDescription(alnmetrics.res.tab)

    pdf(file=file.path(outdir, "mapping-summary.pdf"))
    print(stripplot(TOTAL_READS/1e6 + PF_READS_ALIGNED/1e6  ~ CATEGORY | sampleinfo$description, data=alnmetrics.res.tab, auto.key=list(text=c("Reads", "Aligned")), scales=list(x=list(rot=45)), ylab="Reads (millions)", xlab="Category.", par.settings=simpleTheme(col=c("black","red"), pch=21)))
    dev.off()

    pdf(file=file.path(outdir, "mapping-summary-by-category.pdf"))
    print(stripplot(PCT_PF_READS_ALIGNED ~ sampleinfo$description | CATEGORY,groups=sampleinfo$description, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    if ("sample_prj" %in% names(runinfo)) {
        alnmetrics.res.tab$SAMPLE <- do.call("cbind", strsplit(rownames(alnmetrics.res.tab), "\\."))[1,]
        alnmetrics.res.tab$lane <- do.call("c", runinfo[match(alnmetrics.res.tab$SAMPLE, rownames(runinfo)),]$lane)
        alnmetrics.res.tab$project <- do.call("c", runinfo[match(alnmetrics.res.tab$SAMPLE,rownames(runinfo)),]$sample_prj)
        pdf(file=file.path(outdir, "mapping-summary-by-lane.pdf"))
        print(stripplot(PCT_PF_READS_ALIGNED ~ lane | CATEGORY, groups=project, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45), y=list(rot=45, at=NULL))), par.settings=simpleTheme(pch=19)))
        dev.off()
        pdf(file=file.path(outdir, "mapping-summary-by-project.pdf"))
        print(stripplot(PCT_PF_READS_ALIGNED ~ project | CATEGORY, groups=lane, data=alnmetrics.res.tab, auto.key=list(space="right"), scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
        dev.off()
    }
    alnmetrics.res.tab
}

getInsertmetrics <- function(analysisdir, runinfo, outdir, pattern="*.insert_metrics") {
    sample_prj <- FALSE
    if ("sample_prj" %in% names(runinfo)) {
        sample_prj <- TRUE
        insertmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="insert")
    } else {
        insertmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="insert", sampletable=runinfo[,c("lane", "description")])
    }
    insertmetrics.res.tab <- do.call("rbind", lapply(insertmetrics.res, function(x) {x$metrics}))
    sampleinfo <- getLaneAndDescription(insertmetrics.res.tab, concatenate=sample_prj)
    insertmetrics.res.hist <- as.data.frame(do.call("rbind", lapply(names(insertmetrics.res),
                                                                    function(x) {
                                                                        i.insertsize <- grep("insert_size", colnames(insertmetrics.res[[x]]$histogram))
                                                                        insert_size = data.frame(size=insertmetrics.res[[x]]$histogram[,c(i.insertsize)], frcount=NA, rfcount=NA, tandemcount=NA)
                                                                        i.frcount <- grep("fr_count", colnames(insertmetrics.res[[x]]$histogram))
                                                                        i.rfcount <- grep("rf_count", colnames(insertmetrics.res[[x]]$histogram))
                                                                        i.tandemcount <- grep("tandem_count", colnames(insertmetrics.res[[x]]$histogram))
                                                                        if (length(i.frcount)>0) {
                                                                            insert_size$frcount <- insertmetrics.res[[x]]$histogram[,c(i.frcount)]
                                                                        }

                                                                        if (length(i.rfcount)>0) {
                                                                            insert_size$rfcount <- insertmetrics.res[[x]]$histogram[,c(i.rfcount)]
                                                                        }

                                                                        if (length(i.tandemcount)>0) {
                                                                            insert_size$tandemcount <- insertmetrics.res[[x]]$histogram[,c(i.tandemcount)]
                                                                        }

                                                                        tmp <- cbind(insert_size, sample=x)
                                                                        colnames(tmp) <- c("insert_size", "fr_count", "rf_count", "tandem_count", "sample")
                                                                        tmp
                                                                    })))

    pdf(file=file.path(outdir, "insert-frcount-summary.pdf"))
    print(xyplot(fr_count ~ insert_size | sample, data=insertmetrics.res.hist, xlab="Insert size", ylab="Count", xlim=c(0,1000), main="Insert size distributions", type="l", lwd=2))
    dev.off()
    pdf(file=file.path(outdir, "insert-rfcount-summary.pdf"))
    print(xyplot(rf_count ~ insert_size | sample, data=insertmetrics.res.hist, xlab="Insert size", ylab="Count", xlim=c(0,10000), main="Insert size distributions", type="l", lwd=2))
    dev.off()
    pdf(file=file.path(outdir, "insert-tandem-summary.pdf"))
    print(xyplot(rf_count ~ insert_size | sample, data=insertmetrics.res.hist, xlab="Insert size", ylab="Count", xlim=c(0,10000), main="Insert size distributions", type="l", lwd=2))
    dev.off()
    insertmetrics.res.tab
}

getHsmetrics <- function(analysisdir, runinfo,  outdir, pattern="*.hs_metrics", ...) {
    sample_prj <- FALSE
    if ("sample_prj" %in% names(runinfo)) {
        sample_prj <- TRUE
        hsmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="hs", ...)
    } else {
        hsmetrics.res <- get_hiseq_metrics(analysisdir, runinfo, pattern=pattern, type="hs", sampletable=runinfo[,c("lane", "description")], ...)
    }
    hsmetrics.res.tab <- do.call("rbind", lapply(hsmetrics.res, function(x) {x$metrics}))
    sampleinfo <- getLaneAndDescription(hsmetrics.res.tab, concatenate=sample_prj)

    if ("sample_prj" %in% names(runinfo)) {
        hsmetrics.res.tab$project <- do.call("c", runinfo[match(sampleinfo$description, rownames(runinfo)),]$sample_prj)
    } else  {
        hsmetrics.res.tab$project <- NA
    }

    hsmetrics.res.tab <- cbind(hsmetrics.res.tab, PERCENT_ON_TARGET=hsmetrics.res.tab$FOLD_ENRICHMENT/(hsmetrics.res.tab$GENOME_SIZE / hsmetrics.res.tab$TARGET_TERRITORY)* 100)
    pdf(file=file.path(outdir, "hs-summary.pdf"))
    print(xyplot(FOLD_ENRICHMENT / (GENOME_SIZE/TARGET_TERRITORY) ~ sampleinfo$description, data=hsmetrics.res.tab, xlab="Sample", ylab="Percent on target", main="Percent on target", scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
    dev.off()

    i <- c("ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X")
    hsmetrics.stack <- stack(hsmetrics.res.tab[,i])
    hsmetrics.stack <- cbind(hsmetrics.stack, SAMPLE=sampleinfo$description, project=hsmetrics.res.tab$project)
    hsmetrics.stack$ind <- factor(hsmetrics.stack$ind, levels=levels(hsmetrics.stack$ind)[c(5,3,1,2,4)])
    levels(hsmetrics.stack$ind) <- c("0X", "2X", "10X", "20X", "30X")
    pdf(file=file.path(outdir, "hs-coverage-by-sample.pdf"))
    print(stripplot(values ~ ind | sampleinfo$description, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
    dev.off()

    if ("sample_prj" %in% names(runinfo)) {
        pdf(file=file.path(outdir, "hs-summary-by-project.pdf"))
        print(stripplot(100*FOLD_ENRICHMENT / (GENOME_SIZE/TARGET_TERRITORY) ~ project, data=hsmetrics.res.tab, xlab="Project", ylab="Percent on target", main="Percent on target", scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19)))
        dev.off()

        pdf(file=file.path(outdir, "hs-coverage-by-project.pdf"))
        print(stripplot(values ~ ind | project, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
        dev.off()
        pdf(file=file.path(outdir, "hs-coverage-by-project-bw.pdf"))
        print(bwplot(values ~ ind | project, data=hsmetrics.stack, scales=list(x=list(rot=45)), main="Percentage bases with a given coverage.", xlab="Coverage", ylab="Percentage bases", ylim=c(0,100), par.settings=simpleTheme(pch=19)))
        dev.off()
        if (length(setdiff(rownames(res.df), rownames(hsmetrics.res.tab))) > 0) {
            hsmetrics.res.tab[rownames(res.df[!rownames(res.df) %in% rownames(hsmetrics.res.tab),]),] <- NA
        }
        hsmetrics.res.tab <- hsmetrics.res.tab[match(rownames(res.df), rownames(hsmetrics.res.tab)), ]
    }

    hsmetrics.res.tab
}

## Assumes metrics are named lane_[description|index]
getLaneAndDescription <- function(metrics.tab, concatenate=FALSE) {
    samples <- rownames(metrics.tab)
    lane <- do.call("c", lapply(strsplit(rownames(metrics.tab), "_"), function(x){x[1]}))
    description <- do.call("c", lapply(strsplit(rownames(metrics.tab), "_"), function(x){paste(x[2:length(x)], sep="_", collapse="_")}))
    description <- gsub("\\.[123]$", "", description)
    if (concatenate) {
        description <- paste(lane, description, sep="_")
    }
    data.frame(lane = lane, description = description)
}


