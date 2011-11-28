\name{print.xtable.booktabs}
\alias{print.xtable.booktabs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.xtable.booktabs(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    cmidrule <- ""
    toprule <- "\\toprule "
    if (!is.null(n.cgroup(x))) {
        cmidrule.tmp <- paste(lapply(n.cgroup(x), function(y) {
            paste("\\cmidrule(r){", y[1], "-", y[2], "}")
        }), collapse = " ")
        cmidrule <- paste("\\\\", cmidrule.tmp)
        toprule.tmp <- paste(paste(lapply(n.cgroup(x), function(y) {
            width <- y[2] - y[1] + 1
            paste("\\multicolumn{", width, "}{c}", sep = "")
        }), "{", cgroup(x), "} ", sep = ""), collapse = "& ")
        toprule <- paste(toprule, toprule.tmp, cmidrule)
    }
    print(x, hline.after = NULL, add.to.row = list(pos = list(-1, 
        0, nrow(x)), command = c(toprule, "\\midrule ", "\\bottomrule ")), 
        ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
