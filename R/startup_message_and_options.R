# Package startup message
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Welcome! :)")
}

# Custom options
.onLoad <- function(libname, pkgname) {
    op <- options()  # set up options here
    op.devtools <- list(devtools.path = "~/R-dev", devtools.install.args = "", devtools.name = "Sofia Persson", devtools.desc.author = "Persson Sofia <sofiapersson27@gmail.com>", 
        devtools.desc.license = "What license is it under?", devtools.desc.suggests = NULL, devtools.desc = list())
    toset <- !(names(op.devtools) %in% names(op))
    if (any(toset)) 
        options(op.devtools[toset])
    
    invisible()
}
