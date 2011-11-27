print.xtable.booktabs <-
function(x, ...) {
    cmidrule <- ""
    toprule <- "\\toprule "
    if (!is.null(n.cgroup(x))) {
        cmidrule.tmp <- paste(lapply(n.cgroup(x), function(y) {paste("\\cmidrule(r){", y[1], "-", y[2], "}")}), collapse=" ")
        cmidrule <- paste("\\\\", cmidrule.tmp)
        toprule.tmp <- paste(paste(lapply(n.cgroup(x), function(y) {width <- y[2]-y[1]+1; paste("\\multicolumn{",width, "}{c}", sep="")}), "{", cgroup(x), "} ", sep=""), collapse="& ")
        toprule <- paste(toprule, toprule.tmp, cmidrule)
    }
    print(x,
          hline.after=NULL,
          add.to.row=list(pos=list(-1,0, nrow(x)),
          command=c(
          toprule,
          "\\midrule ",
          "\\bottomrule ")),
          ...)
}

cgroup <-
function(x, ...) UseMethod("cgroup")

`cgroup<-` <-
function(x, value) UseMethod("cgroup<-")

cgroup.xtable <-
function(x, ...) {
 return(attr(x,"cgroup",exact=TRUE))
}

`cgroup<-.xtable` <-
function(x, value) {
    if(!is.null(n.cgroup(x)))
        if(length(n.cgroup(x))!=length(value))
            stop("cgroup must have same length as n.cgroup")
    attr(x, "cgroup") <- value
    return(x)
}
