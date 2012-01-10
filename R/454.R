read454.metrics <- function(indir) {
    infile <- file.path(indir, "454AllControlMetrics.csv")
    lines <- gsub("\"", "", readLines(infile))
    i <- grep("^Key", lines)
    res <- list()
    res$data <- rbind(cbind(region="region1", read.csv(textConnection(lines[(i[1]):(i[1]+2)]))),
                      cbind(region="region2", read.csv(textConnection(lines[(i[2]):(i[2]+2)]))))
    i <- grep("^Accuracy", lines )
    j <- grep("^region", lines )
    ## res$region1$accuracy <- read.csv(textConnection(lines[(i[1]):(j[2]-2)]))
    ## res$region2$accuracy <- read.csv(textConnection(lines[(i[2]):length(lines)]))
    res$accuracy <- rbind(cbind(region="region1", read.csv(textConnection(lines[(i[1]):(j[2]-2)]))),
                          cbind(region="region2", read.csv(textConnection(lines[(i[2]):length(lines)]))))
    return(res)
}
