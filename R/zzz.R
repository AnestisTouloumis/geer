.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("geer", pkgname)
  }

  marginaleffects_classes <- getOption("marginaleffects_model_classes")
  options(
    marginaleffects_model_classes = unique(
      c(marginaleffects_classes, "geer")
    )
  )
}


.onUnload <- function(libpath) {
  marginaleffects_classes <- getOption("marginaleffects_model_classes")
  if (!is.null(marginaleffects_classes)) {
    remaining <- setdiff(marginaleffects_classes, "geer")
    options(
      marginaleffects_model_classes = if (length(remaining)) remaining else NULL
    )
  }

  library.dynam.unload("geer", libpath)
}
