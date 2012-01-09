DEPTHOFCOVERAGE_LABELS = c("cumulative_coverage_counts", "cumulative_coverage_proportions", "interval_statistics", "interval_summary", "statistics", "summary")

read.DepthOfCoverage <- function(prefix, indir="./", ...) {
    X <- lapply(DEPTHOFCOVERAGE_LABELS, function(x) {
        infile <- file.path(indir, paste(prefix, paste("sample_", x, sep=""), sep="."))
        y <- NA
        if (file.exists(infile)) {
            y <- read.table(infile, header=TRUE, fill=TRUE)
        } else {
            infile <- paste(infile, "gz", sep=".")
            y <- read.table(infile, header=TRUE, fill=TRUE)
        }
        if (x == "summary")
            y[y=="N/A"] <- NA
        y
    })
    names(X) <- DEPTHOFCOVERAGE_LABELS
    X
}

setGeneric("plotDepthOfCoverage", function(X, ...) {standardGeneric("plotDepthOfCoverage")})
setMethod("plotDepthOfCoverage",
          signature=signature(
          X="list"
          ),
          function(X, which="cumulative_coverage_counts", pngout=NULL, type="l",  lw=1.3, xlab="coverage", ylab="", legend.pos="topright", relative=FALSE, ...) {
              lab <- match.arg(which, DEPTHOFCOVERAGE_LABELS[1:6])
              Y <- X[[lab]]
              if (relative)
                  Y <- Y/max(Y[,1])
              samples <- rownames(Y)
              n <- dim(Y)[1]
              if (!is.null(pngout))
                  png(pngout)
              matplot(t(Y), main=lab, ylab=ylab, type=type, lw=lw, xlab=xlab, lty=1:n, col=1:n, ...)
              legend(legend.pos, col=1:n, legend=samples, lty=1:n)
              if (!is.null(pngout))
                  dev.off()
          }
          )
setGeneric("plotSampleIntervalSummary", function(X, ...) {standardGeneric("plotSampleIntervalSummary")})
setMethod("plotSampleIntervalSummary",
          signature=signature(
          X="list"
          ),
          function(X, ...) {
          if (exists("X$interval_summary"))
              plotSampleIntervalSummary(X$interval_summary, ...)
          })
setMethod("plotSampleIntervalSummary",
          signature=signature(
          X="data.frame"
          ),
          function(X, which="mean_cvg", plotfn="boxplot", ...) {
              interval_summary <- c("total_cvg","mean_cvg","granular_Q1", "granular_median", "granular_Q3", "%_above_15")
              which <- match.arg(which, interval_summary)
              i <- grep(which, colnames(X))
              boxplot(t(X[,i]), ...)
          })


