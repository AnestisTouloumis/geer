testthat::local_edition(3)

test_that("new_geer and geer construct geer objects from valid lists", {
  stub <- unclass(fit_geewa_pois_indep)

  expect_s3_class(new_geer(stub), "geer")
  expect_s3_class(geer(stub), "geer")
})

test_that("geer rejects invalid inputs", {
  expect_error(geer(42), "'x' must be a list")
  expect_error(geer("a string"), "'x' must be a list")

  incomplete <- list(coefficients = 1, residuals = 1)

  expect_error(
    geer(incomplete),
    "missing required components"
  )
  expect_error(
    geer(incomplete),
    "fitted.values|rank|family"
  )
})


make_valid_geer_input <- function() {
  dat <- data.frame(
    y = c(1, 2, 3),
    x1 = c(0, 1, 0),
    x2 = c(1.0, 0.5, -0.5),
    id = c(1, 1, 2),
    repeated = c(1, 2, 1)
  )

  form <- y ~ x1 + x2
  xmat <- stats::model.matrix(form, dat)
  beta <- c(`(Intercept)` = 0.2, x1 = 0.5, x2 = -0.1)
  fitted <- as.numeric(xmat %*% beta)
  resid <- dat$y - fitted
  V <- diag(c(0.10, 0.20, 0.30))
  dimnames(V) <- list(names(beta), names(beta))

  list(
    coefficients = beta,
    residuals = resid,
    fitted.values = fitted,
    qr = qr(xmat),
    rank = ncol(xmat),
    family = stats::gaussian(),
    linear.predictors = fitted,
    iter = 3L,
    prior.weights = rep(1, nrow(dat)),
    df.residual = nrow(dat) - ncol(xmat),
    y = dat$y,
    x = xmat,
    id = dat$id,
    repeated = dat$repeated,
    converged = TRUE,
    call = quote(geewa(formula = y ~ x1 + x2)),
    formula = form,
    terms = stats::terms(form),
    data = dat,
    offset = NULL,
    control = list(),
    method = "gee",
    naive_covariance = V,
    robust_covariance = V,
    bias_corrected_covariance = V,
    association_structure = "exchangeable",
    alpha = 0.2,
    phi = 1,
    obs_no = nrow(dat),
    clusters_no = length(unique(dat$id)),
    min_cluster_size = 1L,
    max_cluster_size = 2L
  )
}

make_valid_geer_object <- function() {
  new_geer(make_valid_geer_input())
}

test_that("new_geer rejects non-list input", {
  expect_error(new_geer(1), "'x' must be a list")
})

test_that("new_geer rejects missing required components", {
  x <- make_valid_geer_input()
  x$phi <- NULL
  x$phi <- NULL
  x <- x[setdiff(names(x), "phi")]
  expect_error(
    new_geer(x),
    "Input is missing required components: phi"
  )
})

test_that("new_geer adds geer class and fills optional defaults", {
  x <- make_valid_geer_input()
  x <- x[setdiff(names(x), c("na.action", "contrasts", "xlevels"))]

  out <- new_geer(x)

  expect_s3_class(out, "geer")
  expect_null(out$na.action)
  expect_null(out$contrasts)
  expect_equal(out$xlevels, list())
})

test_that("new_geer preserves existing optional fields", {
  x <- make_valid_geer_input()
  x$na.action <- structure(1L, class = "omit")
  x$contrasts <- list(x1 = "contr.treatment")
  x$xlevels <- list(x1 = c("0", "1"))

  out <- new_geer(x)

  expect_identical(out$na.action, x$na.action)
  expect_identical(out$contrasts, x$contrasts)
  expect_identical(out$xlevels, x$xlevels)
})

test_that("validate_geer accepts a valid geer object", {
  obj <- make_valid_geer_object()
  out <- validate_geer(obj)
  expect_identical(out, obj)
})

test_that("validate_geer rejects non-list input", {
  expect_error(validate_geer(1), "a 'geer' object must be a list")
})

test_that("validate_geer rejects objects without geer class", {
  obj <- make_valid_geer_input()
  expect_error(validate_geer(obj), "object must inherit from class 'geer'")
})

test_that("validate_geer checks required fields", {
  obj <- make_valid_geer_object()
  obj$coefficients <- NULL
  obj <- obj[setdiff(names(obj), "coefficients")]

  expect_error(
    validate_geer(obj),
    "object must inherit from class 'geer'"
  )
})

test_that("validate_geer checks coefficients are numeric and named", {
  obj <- make_valid_geer_object()
  obj$coefficients <- as.character(obj$coefficients)
  expect_error(validate_geer(obj), "'coefficients' must be numeric")

  obj <- make_valid_geer_object()
  names(obj$coefficients) <- NULL
  expect_error(validate_geer(obj), "'coefficients' must be a named numeric vector")
})

test_that("validate_geer checks x dimensions and column names", {
  obj <- make_valid_geer_object()
  obj$x <- as.data.frame(obj$x)
  expect_error(validate_geer(obj), "'x' must be a matrix")

  obj <- make_valid_geer_object()
  obj$x <- obj$x[, -1, drop = FALSE]
  expect_error(
    validate_geer(obj),
    "number of columns of 'x' must match length of 'coefficients'"
  )

  obj <- make_valid_geer_object()
  colnames(obj$x)[1] <- "bad"
  expect_error(
    validate_geer(obj),
    "column names of 'x' must match names of 'coefficients'"
  )
})

test_that("validate_geer checks fitted values residuals and observation-aligned components", {
  obj <- make_valid_geer_object()
  obj$fitted.values <- as.character(obj$fitted.values)
  expect_error(validate_geer(obj), "'fitted.values' must be numeric")

  obj <- make_valid_geer_object()
  obj$residuals <- as.character(obj$residuals)
  expect_error(validate_geer(obj), "'residuals' must be numeric")

  obj <- make_valid_geer_object()
  obj$id <- obj$id[-1]
  expect_error(validate_geer(obj), "length of 'id' must match number of observations")

  obj <- make_valid_geer_object()
  obj$y <- obj$y[-1]
  expect_error(validate_geer(obj), "response 'y' must match number of observations")
})

test_that("validate_geer checks obs_no and repeated consistency", {
  obj <- make_valid_geer_object()
  obj$obs_no <- 2.5
  expect_error(validate_geer(obj), "'obs_no' must be a single non-negative integer")

  obj <- make_valid_geer_object()
  obj$obs_no <- 2L
  expect_error(validate_geer(obj), "'obs_no' is inconsistent with observation-level components")

  obj <- make_valid_geer_object()
  obj$repeated <- obj$repeated[-1]
  expect_error(validate_geer(obj), "'repeated' must match number of observations")
})

test_that("validate_geer checks terms call and formula consistency", {
  obj <- make_valid_geer_object()
  obj$terms <- list()
  expect_error(validate_geer(obj), "'terms' must be a terms object")

  obj <- make_valid_geer_object()
  obj$call <- 1
  expect_error(validate_geer(obj), "'call' must be a call")

  obj <- make_valid_geer_object()
  obj$call <- quote(geewa(data = data.frame()))
  expect_error(validate_geer(obj), "'call' must contain a formula component")

  obj <- make_valid_geer_object()
  obj$call <- quote(geewa(formula = y ~ x1))
  expect_error(validate_geer(obj), "'formula' and")
})

test_that("validate_geer checks family structure", {
  obj <- make_valid_geer_object()
  obj$family <- 1
  expect_error(validate_geer(obj), "'family' must be a family object")

  obj <- make_valid_geer_object()
  obj$family$family <- NULL
  expect_error(validate_geer(obj), "'family\\$family' must be a length-1 character value")

  obj <- make_valid_geer_object()
  obj$family$link <- NULL
  expect_error(validate_geer(obj), "'family\\$link' must be a length-1 character value")

  obj <- make_valid_geer_object()
  obj$family$linkinv <- NULL
  expect_error(validate_geer(obj), "'family\\$linkinv' must be a function")
})

test_that("validate_geer checks alpha association structure method and df.residual", {
  obj <- make_valid_geer_object()
  obj$alpha <- "a"
  expect_error(validate_geer(obj), "'alpha' must be numeric")

  obj <- make_valid_geer_object()
  obj$association_structure <- c("a", "b")
  expect_error(validate_geer(obj), "'association_structure' must be a length-1 character value")

  obj <- make_valid_geer_object()
  obj$method <- c("gee", "other")
  expect_error(validate_geer(obj), "'method' must be a length-1 character value")

  obj <- make_valid_geer_object()
  obj$df.residual <- c(1, 2)
  expect_error(validate_geer(obj), "'df.residual' must be a single finite numeric value")
})

test_that("validate_geer checks rank converged and phi", {
  obj <- make_valid_geer_object()
  obj$rank <- -1
  expect_error(validate_geer(obj), "'rank' must be a single non-negative integer")

  obj <- make_valid_geer_object()
  obj$converged <- NA
  expect_error(validate_geer(obj), "'converged' must be a single non-missing logical value")

  obj <- make_valid_geer_object()
  obj$phi <- NA_real_
  expect_error(validate_geer(obj), "'phi' must be a single finite numeric value")
})

test_that("validate_geer checks xlevels and contrasts", {
  obj <- make_valid_geer_object()
  obj$xlevels <- 1
  expect_error(validate_geer(obj), "'xlevels' must be NULL or a named list")

  obj <- make_valid_geer_object()
  obj$xlevels <- list(c("a", "b"))
  expect_error(validate_geer(obj), "'xlevels' must be a named list")

  obj <- make_valid_geer_object()
  obj$contrasts <- 1
  expect_error(validate_geer(obj), "'contrasts' must be NULL or a named list")

  obj <- make_valid_geer_object()
  obj$contrasts <- list("contr.treatment")
  expect_error(validate_geer(obj), "'contrasts' must be a named list")
})

test_that("validate_geer checks cluster counts and cluster-size metadata", {
  obj <- make_valid_geer_object()
  obj$clusters_no <- 0
  expect_error(validate_geer(obj), "'clusters_no' must be a single positive integer")

  obj <- make_valid_geer_object()
  obj$clusters_no <- 3L
  expect_error(validate_geer(obj), "'clusters_no' is inconsistent with 'id'")

  obj <- make_valid_geer_object()
  obj$min_cluster_size <- 0
  expect_error(validate_geer(obj), "'min_cluster_size' must be a single positive integer")

  obj <- make_valid_geer_object()
  obj$max_cluster_size <- 0
  expect_error(validate_geer(obj), "'max_cluster_size' must be a single positive integer")

  obj <- make_valid_geer_object()
  obj$min_cluster_size <- 3L
  obj$max_cluster_size <- 2L
  expect_error(validate_geer(obj), "'min_cluster_size' cannot exceed 'max_cluster_size'")
})

test_that("validate_geer checks covariance matrix structure and dimnames", {
  obj <- make_valid_geer_object()
  obj$naive_covariance <- 1
  expect_error(validate_geer(obj), "'naive_covariance' must be a numeric matrix")

  obj <- make_valid_geer_object()
  obj$robust_covariance <- diag(2)
  expect_error(
    validate_geer(obj),
    "'robust_covariance' must have dimensions compatible with 'coefficients'"
  )

  obj <- make_valid_geer_object()
  V <- obj$bias_corrected_covariance
  rownames(V) <- NULL
  obj$bias_corrected_covariance <- V
  expect_error(
    validate_geer(obj),
    "'bias_corrected_covariance' must have row and column names"
  )

  obj <- make_valid_geer_object()
  V <- obj$robust_covariance
  rownames(V) <- c("a", "b", "c")
  obj$robust_covariance <- V
  expect_error(
    validate_geer(obj),
    "dimnames of 'robust_covariance' must match names of 'coefficients'"
  )
})

test_that("validate_geer checks optional df_adjusted_covariance when present", {
  obj <- make_valid_geer_object()
  V <- diag(length(obj$coefficients))
  dimnames(V) <- list(names(obj$coefficients), names(obj$coefficients))
  obj$df_adjusted_covariance <- V
  expect_identical(validate_geer(obj), obj)

  obj$df_adjusted_covariance <- diag(2)
  expect_error(
    validate_geer(obj),
    "'df_adjusted_covariance' must have dimensions compatible with 'coefficients'"
  )
})
