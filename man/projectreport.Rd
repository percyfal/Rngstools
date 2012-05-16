\name{projectreport}
\alias{projectreport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
projectreport
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
projectreport(analysisdir, run_info, outdir, dupmetrics = list(run = TRUE, pattern = "*.dup_metrics"), insertmetrics = list(run = TRUE, pattern = "*.insert_metrics"), alnmetrics = list(run = TRUE, pattern = "*.align_metrics"), hsmetrics = list(run = TRUE, pattern = "*.hs_metrics"), gcmetrics = list(run = TRUE, pattern = "*.gc_metrics"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{analysisdir}{
%%     ~~Describe \code{analysisdir} here~~
}
  \item{run_info}{
%%     ~~Describe \code{run_info} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
}
  \item{dupmetrics}{
%%     ~~Describe \code{dupmetrics} here~~
}
  \item{insertmetrics}{
%%     ~~Describe \code{insertmetrics} here~~
}
  \item{alnmetrics}{
%%     ~~Describe \code{alnmetrics} here~~
}
  \item{hsmetrics}{
%%     ~~Describe \code{hsmetrics} here~~
}
  \item{gcmetrics}{
%%     ~~Describe \code{gcmetrics} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (analysisdir, run_info, outdir, dupmetrics = list(run = TRUE, 
    pattern = "*.dup_metrics"), insertmetrics = list(run = TRUE, 
    pattern = "*.insert_metrics"), alnmetrics = list(run = TRUE, 
    pattern = "*.align_metrics"), hsmetrics = list(run = TRUE, 
    pattern = "*.hs_metrics"), gcmetrics = list(run = TRUE, pattern = "*.gc_metrics")) 
{
    tmp <- yaml.load_file(run_info)
    if (names(tmp) == "details") {
        tmp <- tmp$details
    }
    runinfo <- as.data.frame(do.call("rbind", lapply(tmp, function(x) {
        y <- as.data.frame(rbind(x))
        rownames(y) <- paste(x$lane, x$description, sep = "_")
        y
    })))
    fcinfo <- do.call("rbind", lapply(runinfo$files, function(x) {
        m <- regexec("(^[0-9]+).*_([0-9]{6})_([0-9A-Z]+)_.*", 
            basename(x)[[1]])
        y <- regmatches(basename(x)[[1]], m)
        y[[1]][2:4]
    }))
    colnames(fcinfo) <- c("fc_lane", "date", "flowcell")
    runinfo <- cbind(runinfo, fcinfo)
    samples <- runinfo$description
    res.list <- list(dupmetrics = NA, insertmetrics = NA, alnmetrics = NA, 
        alnmetrics_pair = NA, alnmetrics_first_of_pair = NA, 
        alnmetrics_second_of_pair = NA)
    res.df <- NULL
    if (dupmetrics$run) {
        cat("Adding dupmetrics...\n")
        dupmetrics.res.tab <- getDupmetrics(analysisdir, runinfo, 
            outdir, dupmetrics$pattern)
        res.list$dupmetrics <- dupmetrics.res.tab
        res.df <- dupmetrics.res.tab
        reportfile <- file.path(outdir, "dupmetrics.txt")
        write.table(file = reportfile, dupmetrics.res.tab, sep = "\t", 
            row.names = FALSE)
    }
    if (alnmetrics$run) {
        cat("Adding alnmetrics...\n")
        alnmetrics.res.tab <- getAlnmetrics(analysisdir, runinfo, 
            outdir, pattern = alnmetrics$pattern)
        res.list$alnmetrics <- alnmetrics.res.tab
        res.list$alnmetrics_pair = alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY == 
            "PAIR", ]
        colnames(res.list$alnmetrics_pair) <- paste("aln", colnames(res.list$alnmetrics_pair), 
            "3", sep = "_")
        res.list$alnmetrics_first_of_pair = alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY == 
            "FIRST_OF_PAIR", ]
        colnames(res.list$alnmetrics_first_of_pair) <- paste("aln", 
            colnames(res.list$alnmetrics_first_of_pair), "1", 
            sep = "_")
        res.list$alnmetrics_second_of_pair = alnmetrics.res.tab[alnmetrics.res.tab$CATEGORY == 
            "SECOND_OF_PAIR", ]
        colnames(res.list$alnmetrics_second_of_pair) <- paste("aln", 
            colnames(res.list$alnmetrics_second_of_pair), "2", 
            sep = "_")
        res.df <- cbind(res.df, res.list$alnmetrics_pair, res.list$alnmetrics_first_of_pair, 
            res.list$alnmetrics_second_of_pair)
        reportfile <- file.path(outdir, "alnmetrics.txt")
        write.table(file = reportfile, alnmetrics.res.tab, sep = "\t", 
            row.names = FALSE)
    }
    if (insertmetrics$run) {
        cat("Adding insertmetrics...\n")
        insertmetrics.res.tab <- getInsertmetrics(analysisdir, 
            runinfo, outdir, pattern = insertmetrics$pattern)
        res.list$insertmetrics = insertmetrics.res.tab
        res.df <- cbind(res.df, insertmetrics.res.tab)
        reportfile <- file.path(outdir, "insertmetrics.txt")
        write.table(file = reportfile, insertmetrics.res.tab, 
            sep = "\t", row.names = FALSE)
    }
    if (hsmetrics$run) {
        cat("Adding hsmetrics...\n")
        hsmetrics.res.tab <- getHsmetrics(analysisdir, runinfo, 
            outdir, pattern = hsmetrics$pattern)
        res.list$hsmetrics = hsmetrics.res.tab
        res.df <- cbind(res.df, hsmetrics.res.tab)
        reportfile <- file.path(outdir, "hsmetrics.txt")
        write.table(file = reportfile, hsmetrics.res.tab, sep = "\t", 
            row.names = FALSE)
    }
    res.df <- cbind(res.df, runinfo[match(rownames(res.df), rownames(runinfo)), 
        ])
    res.df$lane <- as.integer(res.df$lane)
    list(res.list, res.df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
