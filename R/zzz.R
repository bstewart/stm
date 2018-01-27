.onAttach <- function(libname, pkgname) {
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pd$Version, " (",
                    pd$Date, ") successfully loaded. See ?", pkgname, " for help. \n Papers, resources, and other materials at structuraltopicmodel.com");
}