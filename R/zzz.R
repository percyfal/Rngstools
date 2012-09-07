## zzz.R (2011-09-05)

##   Library Loading

.RngstoolsEnv <- new.env()

.onLoad <- function(libname, pkgname) {
    .setRngstoolsOptions(pkgname)
}
