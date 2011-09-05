\name{read.diBayes.gff}
\alias{read.diBayes.gff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  read.diBayes.gff
}
\description{
  Read a diBayes gff output file
}
\usage{
read.diBayes.gff(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function(file) {
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
    gff$titv <- mapply(function(x,y) {.titv.df[x,y]}, gff$genotype, gff$reference)
    gff
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
