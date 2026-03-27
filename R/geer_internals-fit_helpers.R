build_geer_model_frame <- function(model_call, env = parent.frame()) {
  model_call$drop.unused.levels <- TRUE
  keep <- c("formula", "data", "id", "repeated", "weights", "offset")
  model_call <- model_call[c(1L, match(keep, names(model_call), 0L))]
  model_call[[1L]] <- quote(model.frame)

  enclos <- env
  if (!is.null(model_call$formula)) {
    formula_obj <- tryCatch(eval(model_call$formula, envir = env), error = function(e) NULL)
    if (inherits(formula_obj, "formula") && !is.null(environment(formula_obj))) {
      enclos <- environment(formula_obj)
    }
  }

  eval(model_call, envir = env, enclos = enclos)
}


extract_geer_offset <- function(model_frame, y_length) {
  offset <- model.offset(model_frame)
  if (is.null(offset)) {
    offset <- rep.int(0, y_length)
  } else {
    if (!is.numeric(offset)) {
      stop("'offset' should be a numeric vector", call. = FALSE)
    }
    if (length(offset) == 1L) {
      offset <- rep.int(as.numeric(offset), y_length)
    }
    if (length(offset) != y_length) {
      stop("response variable and 'offset' are not of same length", call. = FALSE)
    }
    offset <- as.double(offset)
  }
  offset
}


build_geer_design_matrix <- function(model_frame) {
  model_terms <- attr(model_frame, "terms")
  model_matrix <- model.matrix(model_terms, model_frame)
  xnames <- colnames(model_matrix)
  if (length(xnames) == 1L) {
    model_matrix <- matrix(model_matrix, ncol = 1L)
  }
  qr_model_matrix <- qr(model_matrix)
  if (qr_model_matrix$rank < ncol(model_matrix)) {
    stop("rank-deficient model matrix", call. = FALSE)
  }
  list(
    terms = model_terms,
    x = model_matrix,
    xnames = xnames,
    qr = qr_model_matrix,
    assign = attr(model_matrix, "assign"),
    contrasts = attr(model_matrix, "contrasts")
  )
}


normalize_geer_control <- function(control) {
  if (is.null(control)) {
    return(geer_control())
  }
  if (is.list(control) &&
      !is.null(control$maxiter) &&
      !is.null(control$tolerance)) {
    return(control)
  }
  do.call("geer_control", control)
}
