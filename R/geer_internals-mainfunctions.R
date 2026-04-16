valid_methods <- c(
  "gee",
  "brgee-naive", "brgee-robust", "brgee-empirical",
  "bcgee-naive", "bcgee-robust", "bcgee-empirical",
  "pgee-jeffreys", "opgee-jeffreys", "hpgee-jeffreys"
)
valid_corstrs <- c(
  "independence", "exchangeable", "ar1",
  "m-dependent", "unstructured", "fixed"
)
valid_orstrs <- c(
  "independence", "exchangeable", "unstructured", "fixed"
)
valid_families <- c(
  "gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian",
  "quasi", "quasibinomial", "quasipoisson"
)
valid_links <- c(
  "logit", "probit", "cauchit", "cloglog", "identity", "log",
  "sqrt", "1/mu^2", "inverse"
)


normalize_family <- function(family) {
  if (is.character(family)) {
    family_fun <- tryCatch(
      match.fun(family),
      error = function(e) NULL
    )
    if (is.null(family_fun)) {
      stop(
        "'family' must be a family object, a family function, or the name of one",
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
  check_choice(family$family, valid_families, "family")
  link <- as.character(family$link)
  check_choice(link, valid_links, "link")
  family
}

extract_geer_response_weights <- function(model_frame, family) {
  is_binomial <- identical(family$family, "binomial")
  y <- model.response(model_frame, "any")
  if (is.null(y)) stop("response variable not found", call. = FALSE)
  y_raw <- y
  if (is_binomial) {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop(
          "for binomial models, a factor response must have exactly two levels",
          call. = FALSE
        )
      }
      y <- as.numeric(y == levels(y)[2L])
    } else if (is.character(y)) {
      yfac <- factor(y)
      if (nlevels(yfac) != 2L) {
        stop(
          "for binomial models, a character response must have exactly two distinct values",
          call. = FALSE
        )
      }
      y <- as.numeric(yfac == levels(yfac)[2L])
    } else if (is.matrix(y) && ncol(y) != 2L) {
      stop(
        "for binomial models, a matrix response must have exactly two columns",
        call. = FALSE
      )
    }
  }
  weights <- as.vector(model.weights(model_frame))
  if (is.null(weights)) {
    weights <- rep.int(1, length(y))
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
    if (length(weights) != length(y)) {
      stop("'weights' and the response must have the same length", call. = FALSE)
    }
  }
  weights <- as.numeric(weights)

  if (is_binomial && is.matrix(y) && ncol(y) == 2L) {
    trials <- rowSums(y)
    if (any(!is.finite(trials)) || any(trials <= 0)) {
      stop(
        "for binomial matrix responses, row sums (trials) must be positive and finite",
        call. = FALSE
      )
    }
    weights <- weights * trials
    y <- y[, 1L] / trials
  }
  y <- as.numeric(y)
  if (any(!is.finite(y))) {
    stop("response variable contains non-finite values", call. = FALSE)
  }
  if (is_binomial && !is.matrix(y_raw)) {
    if (any(y < 0 | y > 1)) {
      stop(
        "for binomial models, the response must be coded as 0/1, proportions in [0, 1], or a two-column matrix",
        call. = FALSE
      )
    }
  }
  list(y = y, weights = weights)
}


extract_geer_id_repeated <- function(model_frame, y_length) {
  id_raw <- model.extract(model_frame, "id")
  if (is.null(id_raw)) stop("'id' not found", call. = FALSE)
  if (anyNA(id_raw)) stop("'id' cannot contain missing values", call. = FALSE)
  id <- as.numeric(factor(id_raw))
  if (length(id) != y_length) {
    stop("response variable and 'id' are not of same length", call. = FALSE)
  }
  repeated <- model.extract(model_frame, "repeated")
  if (is.null(repeated)) {
    repeated <- ave(id, id, FUN = seq_along)
  } else {
    if (anyNA(repeated)) stop("'repeated' cannot contain missing values", call. = FALSE)
    repeated <- as.numeric(factor(repeated))
  }
  if (length(repeated) != y_length) {
    stop("response variable and 'repeated' are not of same length", call. = FALSE)
  }
  if (any(unlist(lapply(split(repeated, id), duplicated)))) {
    stop("'repeated' does not have unique values per 'id'", call. = FALSE)
  }
  list(id = id, repeated = repeated)
}



normalize_phi <- function(phi_fixed, phi_value) {
  if (!is.logical(phi_fixed) || length(phi_fixed) != 1L || is.na(phi_fixed)) {
    stop("'phi_fixed' must be a single non-missing logical value", call. = FALSE)
  }
  if (phi_fixed) {
    phi_value <- as.numeric(phi_value)
    if (length(phi_value) != 1L || !is.finite(phi_value) || phi_value <= 0) {
      stop(
        "'phi_value' must be a single positive number when 'phi_fixed = TRUE'",
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
    stop("'use_p' must be a single non-missing logical value", call. = FALSE)
  }
  use_p
}

