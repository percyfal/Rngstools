\name{read.DepthOfCoverage}
\alias{read.DepthOfCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  read.DepthOfCoverage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.DepthOfCoverage(prefix, indir = "./", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
}
  \item{indir}{
%%     ~~Describe \code{indir} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (prefix, indir = "./", ...) 
{
    labels <- c("sample_cumulative_coverage_counts", "sample_cumulative_coverage_proportions", 
        "sample_interval_statistics", "sample_interval_summary", 
        "sample_statistics", "sample_summary")
    X <- lapply(labels, function(x) {
        infile <- file.path(indir, paste(prefix, x, sep = "."))
        y <- NA
        if (file.exists(infile)) {
            y <- read.table(infile, header = TRUE, fill = TRUE)
        }
        else {
            infile <- paste(infile, "gz", sep = ".")
            y <- read.table(infile, header = TRUE, fill = TRUE)
        }
        y
    })
    names(X) <- labels
    X
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
