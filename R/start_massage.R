##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("tinyarray v 2.4.2  welcome to use tinyarray!
If you use tinyarray in published research, please acknowledgements:
We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee,especially Xiaojie Sun, for generously sharing their experience and codes.
")
}
