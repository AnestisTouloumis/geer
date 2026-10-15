build_geer_model_frame <- function(model_call, env = parent.frame()) {
  model_call$drop.unused.levels <- TRUE
  keep <- c(
    "formula", "data", "subset", "na.action",
    "id", "repeated", "weights", "offset"
  )
  model_call <- model_call[c(1L, match(keep, names(model_call), 0L))]
  model_call[[1L]] <- quote(stats::model.frame)
  enclos <- env
  if (!is.null(model_call$formula)) {
    formula_obj <- tryCatch(
      eval(model_call$formula, envir = env),
      error = function(e) NULL
    )
    if (inherits(formula_obj, "formula") &&
        !is.null(environment(formula_obj))) {
      enclos <- environment(formula_obj)
    }
  }
  model_call <- recycle_geer_scalar_offset(model_call, env, enclos)
  eval(model_call, envir = env, enclos = enclos)
}


recycle_geer_scalar_offset <- function(model_call, env, enclos) {
  if (is.null(model_call$offset)) {
    return(model_call)
  }
  data_arg <- if (is.null(model_call$data)) {
    NULL
  } else {
    eval(model_call$data, envir = env, enclos = enclos)
  }
  offset_value <- tryCatch(
    eval(model_call$offset, envir = data_arg, enclos = enclos),
    error = function(e) NULL
  )
  if (!is.numeric(offset_value) || length(offset_value) != 1L) {
    return(model_call)
  }
  frame_call <- model_call[
    c(1L, match(c("formula", "data"), names(model_call), 0L))
  ]
  frame_call$na.action <- quote(stats::na.pass)
  n_rows <- NROW(eval(frame_call, envir = env, enclos = enclos))
  model_call$offset <- rep.int(as.numeric(offset_value), n_rows)
  model_call
}


extract_geer_offset <- function(model_frame, y_length) {
  offset <- stats::model.offset(model_frame)
  if (is.null(offset)) {
    offset <- rep.int(0, y_length)
  } else {
    if (!is.numeric(offset)) {
      stop("'offset' must be a numeric vector", call. = FALSE)
    }
    if (length(offset) == 1L) {
      offset <- rep.int(as.numeric(offset), y_length)
    }
    if (length(offset) != y_length) {
      stop(
        "'offset' must have one value per observation",
        call. = FALSE
      )
    }
    if (anyNA(offset) || any(!is.finite(offset))) {
      stop("'offset' must be finite", call. = FALSE)
    }
    offset <- as.double(offset)
  }
  offset
}


build_geer_design_matrix <- function(model_frame) {
  model_terms <- attr(model_frame, "terms")
  model_matrix <- stats::model.matrix(model_terms, model_frame)
  xnames <- colnames(model_matrix)
  if (ncol(model_matrix) == 0L) {
    stop(
      paste0(
        "the model matrix has no columns; include at least one term or an ",
        "intercept"
      ),
      call. = FALSE
    )
  }
  qr_model_matrix <- qr(model_matrix)
  if (qr_model_matrix$rank < ncol(model_matrix)) {
    aliased <- xnames[qr_model_matrix$pivot[
      seq.int(qr_model_matrix$rank + 1L, ncol(model_matrix))
    ]]
    stop(
      "rank-deficient model matrix: the column(s) ",
      paste0("'", aliased, "'", collapse = ", "),
      " are linearly dependent on the others; remove them from the formula",
      call. = FALSE
    )
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
  if (!is.list(control)) {
    stop(
      "'control' must be NULL or a list, typically created by geer_control()",
      call. = FALSE
    )
  }
  control <- utils::modifyList(geer_control(), control)
  do.call(geer_control, control)
}


normalize_family <- function(family) {
  if (is.character(family)) {
    family_fun <- tryCatch(
      match.fun(family),
      error = function(e) NULL
    )
    if (is.null(family_fun)) {
      stop(
        paste0(
          "'family' must be a family object, a family function, or the name ",
          "of one"
        ),
        call. = FALSE
      )
    }
    family <- family_fun()
  } else if (is.function(family)) {
    family <- family()
  }
  if (!is.list(family) || is.null(family$family) || is.null(family$link)) {
    stop("'family' must be a valid family object", call. = FALSE)
  }
  check_choice(family$family, geer_family_choices, "family")
  link <- as.character(family$link)
  check_choice(link, geer_link_choices, "link")
  family
}


check_geer_response_support <- function(y, family) {
  family_name <- family$family
  if (identical(family_name, "quasi") && !is.null(family$varfun)) {
    if (family$varfun %in% c("mu", "mu^2", "mu^3") && any(y < 0)) {
      stop(
        sprintf(
          "negative values not allowed for quasi(variance = \"%s\")",
          family$varfun
        ),
        call. = FALSE
      )
    }
    return(invisible(NULL))
  }
  message <- switch(
    family_name,
    poisson = if (any(y < 0)) {
      "negative values not allowed for the 'Poisson' family"
    },
    quasipoisson = if (any(y < 0)) {
      "negative values not allowed for the 'quasiPoisson' family"
    },
    Gamma = if (any(y <= 0)) {
      "non-positive values not allowed for the 'Gamma' family"
    },
    inverse.gaussian = if (any(y <= 0)) {
      "positive values only are allowed for the 'inverse.gaussian' family"
    },
    NULL
  )
  if (!is.null(message)) {
    stop(message, call. = FALSE)
  }
  invisible(NULL)
}


extract_geer_response_weights <- function(model_frame, family) {
  is_binomial_like <-
    identical(family$family, "binomial") ||
    identical(family$family, "quasibinomial") ||
    (identical(family$family, "quasi") &&
     !is.null(family$varfun) &&
     identical(family$varfun, "mu(1-mu)"))
  y <- stats::model.response(model_frame, "any")
  if (is.null(y)) {
    stop("the formula must have a response", call. = FALSE)
  }
  if (is_binomial_like) {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop(
          paste0(
            "for binomial-type models, a factor response must have exactly ",
            "two levels"
          ),
          call. = FALSE
        )
      }
      y <- as.numeric(y == levels(y)[2L])
    } else if (is.character(y)) {
      yfac <- factor(y)
      if (nlevels(yfac) != 2L) {
        stop(
          paste0(
            "for binomial-type models, a character response must have exactly ",
            "two distinct values"
          ),
          call. = FALSE
        )
      }
      y <- as.numeric(yfac == levels(yfac)[2L])
    } else if (is.matrix(y) && ncol(y) != 2L) {
      stop(
        paste0(
          "for binomial-type models, a matrix response must have exactly two ",
          "columns"
        ),
        call. = FALSE
      )
    }
  }
  if (!is_binomial_like) {
    if (is.factor(y) || is.character(y)) {
      stop(
        "the response must be numeric unless the family is binomial-type",
        call. = FALSE
      )
    }
    if (is.matrix(y) && ncol(y) > 1L) {
      stop(
        paste0(
          "a matrix response must have one column unless the family is ",
          "binomial-type"
        ),
        call. = FALSE
      )
    }
  }
  n_obs <- if (is.matrix(y)) nrow(y) else length(y)
  weights <- as.vector(stats::model.weights(model_frame))
  if (is.null(weights)) {
    weights <- rep.int(1, n_obs)
  } else {
    if (!is.numeric(weights)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    }
    if (anyNA(weights) || any(!is.finite(weights))) {
      stop("'weights' must be finite", call. = FALSE)
    }
    if (any(weights <= 0)) {
      stop("'weights' must be strictly positive", call. = FALSE)
    }
    if (length(weights) != n_obs) {
      stop(
        "'weights' must have one value per observation",
        call. = FALSE
      )
    }
  }
  weights <- as.numeric(weights)
  if (is_binomial_like && is.matrix(y) && ncol(y) == 2L) {
    if (anyNA(y) || any(!is.finite(y))) {
      stop(
        "for binomial matrix responses, all entries must be finite",
        call. = FALSE
      )
    }
    if (any(y < 0)) {
      stop(
        "for binomial matrix responses, all counts must be non-negative",
        call. = FALSE
      )
    }
    trials <- rowSums(y)
    if (any(trials <= 0)) {
      stop(
        "for binomial matrix responses, row sums (trials) must be positive",
        call. = FALSE
      )
    }
    weights <- weights * trials
    y <- y[, 1L] / trials
  }
  y <- as.numeric(y)
  if (any(!is.finite(y))) {
    stop("the response must contain only finite values", call. = FALSE)
  }
  if (is_binomial_like) {
    if (any(y < 0 | y > 1)) {
      stop(
        paste0(
          "for binomial-type models, the response must be coded as 0/1, ",
          "proportions in [0, 1], or a two-column matrix"
        ),
        call. = FALSE
      )
    }
  } else {
    check_geer_response_support(y, family)
  }
  list(y = y, weights = weights)
}


extract_geer_id_repeated <- function(model_frame,
                                     y_length,
                                     warn_character = TRUE) {
  id_raw <- stats::model.extract(model_frame, "id")
  if (is.null(id_raw)) {
    stop("'id' must be supplied", call. = FALSE)
  }
  if (anyNA(id_raw)) {
    stop("'id' must not contain missing values", call. = FALSE)
  }
  id <- as.numeric(factor(id_raw))
  if (length(id) != y_length) {
    stop("'id' must have one value per observation", call. = FALSE)
  }
  repeated <- stats::model.extract(model_frame, "repeated")
  if (is.null(repeated)) {
    repeated <- stats::ave(id, id, FUN = seq_along)
  } else {
    if (anyNA(repeated)) {
      stop("'repeated' must not contain missing values", call. = FALSE)
    }
    if (warn_character && is.character(repeated)) {
      warning(
        paste0(
          "'repeated' is a character vector and is sorted alphabetically ",
          "(\"10\" before \"2\"); use a numeric or factor variable to set the ",
          "order"
        ),
        call. = FALSE
      )
    }
    repeated <- as.numeric(factor(repeated))
  }
  if (length(repeated) != y_length) {
    stop(
      "'repeated' must have one value per observation",
      call. = FALSE
    )
  }
  if (any(unlist(lapply(split(repeated, id), duplicated)))) {
    stop(
      "'repeated' must have distinct values within each cluster",
      call. = FALSE
    )
  }
  list(id = id, repeated = repeated)
}


normalize_phi <- function(phi_fixed, phi_value) {
  if (!is.logical(phi_fixed) || length(phi_fixed) != 1L || is.na(phi_fixed)) {
    stop(
      "'phi_fixed' must be TRUE or FALSE",
      call. = FALSE
    )
  }
  if (phi_fixed) {
    phi_value <- as.numeric(phi_value)
    if (length(phi_value) != 1L || !is.finite(phi_value) || phi_value <= 0) {
      stop(
        "'phi_value' must be a positive finite number when 'phi_fixed = TRUE'",
        call. = FALSE
      )
    }
  } else {
    phi_value <- 1
  }
  list(phi_fixed = phi_fixed, phi_value = phi_value)
}


normalize_use_p <- function(use_p) {
  if (!is.logical(use_p) || length(use_p) != 1L || is.na(use_p)) {
    stop("'use_p' must be TRUE or FALSE", call. = FALSE)
  }
  use_p
}


prepare_geer_inputs <- function(mcall, family, env, control, method) {
  model_frame <- build_geer_model_frame(mcall, env = env)
  response_weights <- extract_geer_response_weights(model_frame, family)
  y <- response_weights$y
  weights <- response_weights$weights
  id_repeated <- extract_geer_id_repeated(model_frame, length(y))
  id <- id_repeated$id
  repeated <- id_repeated$repeated
  offset <- extract_geer_offset(model_frame, y_length = length(y))
  design <- build_geer_design_matrix(model_frame)
  model_matrix <- design$x
  ord <- order(id, repeated)
  y <- y[ord]
  model_matrix <- model_matrix[ord, , drop = FALSE]
  weights <- weights[ord]
  offset <- offset[ord]
  id <- id[ord]
  repeated <- repeated[ord]
  attr(model_matrix, "assign") <- design$assign
  attr(model_matrix, "contrasts") <- design$contrasts
  control <- normalize_geer_control(control)
  method <- as.character(method)
  check_choice(method, geer_method_choices, "method")
  list(
    row_order = ord,
    model_frame = model_frame,
    y = y,
    weights = weights,
    id = id,
    repeated = repeated,
    offset = offset,
    model_matrix = model_matrix,
    model_terms = design$terms,
    xnames = design$xnames,
    qr_model_matrix = design$qr,
    control = control,
    method = method
  )
}
