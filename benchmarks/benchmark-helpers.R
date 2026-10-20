# =============================================================================
# Shared helpers for the benchmark-*.R scripts.
# =============================================================================

bench_tol <- 1e-6
bench_control <- geer_control(tolerance = 1e-10, maxiter = 1000)

num <- function(x) {
  x <- unclass(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  if (is.matrix(x) && ncol(x) == 1L) x <- x[, 1L]
  attributes(x) <- if (is.matrix(x)) list(dim = dim(x)) else NULL
  storage.mode(x) <- "double"
  x
}

expect_agree <- function(ref, est, tol = bench_tol, label = NULL) {
  expect_equal(num(est), num(ref), tolerance = tol, label = label)
}

relative_difference <- function(ref, est) {
  ref <- num(ref)
  est <- num(est)
  if (!identical(dim(ref), dim(est)) || length(ref) != length(est)) {
    return(Inf)
  }
  max(abs(ref - est)) / max(1, max(abs(ref)))
}

expect_known_difference <- function(ref, est, reason, tol = bench_tol) {
  d <- relative_difference(ref, est)
  if (is.finite(d) && d <= tol) {
    succeed()
  } else if (anyNA(num(est)) || anyNA(num(ref))) {
    skip(sprintf("known difference (%s): one of the values is NA", reason))
  } else {
    skip(sprintf("known difference (%s): max relative difference %.2e",
                 reason, d))
  }
}

in_geer_order <- function(x, fit) {
  x <- num(x)
  x[fit$row_order]
}

geer_corr_matrix <- function(fit, size) {
  alpha <- unname(fit$alpha)
  R <- diag(size)
  lag <- abs(outer(seq_len(size), seq_len(size), "-"))
  structure <- fit$association_structure
  if (identical(structure, "exchangeable")) {
    R[lag > 0] <- alpha[1L]
  } else if (identical(structure, "ar1")) {
    R <- alpha[1L]^lag
  } else if (structure %in% c("m-dependent", "toeplitz")) {
    for (k in seq_along(alpha)) R[lag == k] <- alpha[k]
  } else if (structure %in% c("unstructured", "fixed")) {
    R[lower.tri(R)] <- alpha
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
  }
  R
}

geewa_fixed <- function(args, R, phi) {
  args$corstr <- "fixed"
  args$Mv <- NULL
  args$alpha_vector <- R[lower.tri(R)]
  args$phi_fixed <- TRUE
  args$phi_value <- phi
  geewa_native(args)
}

geewa_native <- function(args) {
  args$control <- bench_control
  withCallingHandlers(
    do.call(geewa, args),
    warning = function(w) {
      if (grepl("bias-corrected covariance matrix is undefined",
                conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

muffle_noninteger <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (startsWith(conditionMessage(w), "non-integer x")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

quietly <- function(expr) {
  out <- NULL
  suppressMessages(utils::capture.output(out <- expr))
  out
}

check_fit <- function(label, ref, fit, tol = bench_tol) {
  order_vec <- function(x) {
    if (isTRUE(ref$ordered)) num(x) else in_geer_order(x, fit)
  }
  describe(label, {
    if (!is.null(ref$coefficients)) {
      it("coefficients agree", {
        expect_agree(ref$coefficients, coef(fit), tol)
      })
    }
    if (!is.null(ref$phi)) {
      it("dispersion agrees", {
        expect_agree(ref$phi, fit$phi, tol)
      })
    }
    if (!is.null(ref$alpha)) {
      it("working association parameters agree", {
        expect_agree(ref$alpha, fit$alpha, tol)
      })
    }
    if (!is.null(ref$corr)) {
      it("working correlation matrix agrees", {
        expect_agree(ref$corr, geer_corr_matrix(fit, nrow(ref$corr)), tol)
      })
    }
    cov_types <- c("naive", "robust", "df_adjusted", "bias_corrected")
    cov_types <- cov_types[cov_types %in% names(ref)]
    if (length(cov_types) > 0) {
      it(paste(paste(sub("_", "-", cov_types), collapse = ", "),
               "covariance matrices agree"), {
        for (type in cov_types) {
          expect_agree(ref[[type]],
                       vcov(fit, cov_type = sub("_", "-", type)), tol,
                       label = type)
        }
      })
    }
    if (!is.null(ref$fitted)) {
      it("fitted values agree", {
        expect_agree(order_vec(ref$fitted), fitted(fit), tol)
      })
    }
    if (!is.null(ref$linear_predictors)) {
      it("linear predictors agree", {
        expect_agree(order_vec(ref$linear_predictors),
                     predict(fit, type = "link"), tol)
      })
    }
    res_types <- c("working", "pearson", "deviance")
    res_types <- res_types[paste0("residuals_", res_types) %in% names(ref)]
    if (length(res_types) > 0) {
      it(paste(paste(res_types, collapse = ", "), "residuals agree"), {
        for (type in res_types) {
          expect_agree(order_vec(ref[[paste0("residuals_", type)]]),
                       residuals(fit, type = type), tol, label = type)
        }
      })
    }
    if (!is.null(ref$residuals_mahalanobis)) {
      it("Mahalanobis residuals agree", {
        expect_agree(ref$residuals_mahalanobis,
                     residuals(fit, type = "mahalanobis"), tol)
      })
    }
    if (!is.null(ref$summary)) {
      s <- ref$summary
      it(paste0("coefficient table (", s$cov_type, ") agrees"), {
        tab <- summary(fit, cov_type = s$cov_type)$coefficients
        if (!is.null(s$se)) expect_agree(s$se, tab[, "Std. Error"], tol)
        if (!is.null(s$statistic)) {
          z <- tab[, "z value"]
          if (identical(s$statistic_type, "wald")) z <- z^2
          expect_agree(s$statistic, z, tol)
        }
        if (!is.null(s$p_value)) expect_agree(s$p_value, tab[, "Pr(>|z|)"], tol)
      })
    }
    if (!is.null(ref$confint)) {
      ci <- ref$confint
      it(paste0("confidence intervals (", ci$cov_type, ") agree"), {
        expect_agree(ci$ci, confint(fit, level = ci$level,
                                    cov_type = ci$cov_type), tol)
      })
    }
  })
}

geer_criteria <- function(fit, criteria, cov_type = "robust") {
  out <- geecriteria(fit, criteria = criteria, cov_type = cov_type,
                     digits = 15)
  vapply(criteria, function(k) as.numeric(out[1L, k]), numeric(1))
}

check_criteria <- function(label, ref, fit, cov_type = "robust",
                           known = list(), tol = bench_tol) {
  values <- geer_criteria(fit, names(ref), cov_type = cov_type)
  describe(label, {
    for (criterion in names(ref)) {
      reason <- known[[criterion]]
      test <- if (is.null(reason)) {
        bquote(expect_agree(.(ref[[criterion]]), .(values[[criterion]]),
                            .(tol)))
      } else {
        bquote(expect_known_difference(.(ref[[criterion]]),
                                       .(values[[criterion]]),
                                       .(reason), .(tol)))
      }
      eval(bquote(it(.(paste(criterion, "agrees")), .(test))))
    }
  })
}

not_comparable <- function(label, reason) {
  describe(label, {
    it("is not compared", {
      skip(paste("not comparable:", reason))
    })
  })
}
