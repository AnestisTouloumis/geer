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
