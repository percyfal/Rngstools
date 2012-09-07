.setRngstoolsOptions <- function(pkgname) {
    RngstoolsOpt <- list(dec=".")
    class(RngstoolsOpt) <- "RngstoolsOptions"
    options("Rngstools"=RngstoolsOpt)

}
