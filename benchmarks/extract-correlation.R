# =============================================================================
# extract_corr_matrix()
#
# Reconstructs the full working correlation matrix from a fitted GEE object.
# Supported classes: geeglm (geepack), gee (gee), geem2 (mmmgee),
#                    glmgee (glmtoolbox).
#
# Arguments:
#   fit      — a fitted GEE object
#   max_size — integer; dimension of the desired square matrix.
#              Required for geeglm, gee, and geem2 where no full matrix is
#              stored directly. For glmgee, inferred automatically from
#              fit$corr.
#
# Returns:
#   A symmetric max_size x max_size correlation matrix with ones on the
#   diagonal, suitable for:
#     geewa(..., corstr = "fixed", alpha_vector = R[lower.tri(R)])
#
# geepack slot reference (from geese.R source):
#   fit$geese$alpha    — compact correlation parameters
#   fit$geese$gamma    — dispersion phi
#   fit$corstr         — correlation structure string
#   fit$geese$vbeta    — robust vcov
#   fit$geese$vbeta.naiv — naive vcov (unscaled; multiply by gamma for phi-scaled)
# =============================================================================

extract_corr_matrix <- function(fit, max_size) {

  # ---------------------------------------------------------------------------
  # glmtoolbox::glmgee — full matrix stored directly in fit$corr
  # ---------------------------------------------------------------------------
  if (inherits(fit, "glmgee")) {
    R <- fit$corr
    if (!is.matrix(R))
      stop("Could not extract correlation matrix from glmgee object.",
           call. = FALSE)
    return(R)
  }

  # ---------------------------------------------------------------------------
  # gee::gee — full matrix stored directly in fit$working.correlation
  # ---------------------------------------------------------------------------
  if (inherits(fit, "gee") && !inherits(fit, "geeglm")) {
    R <- fit$working.correlation
    if (!is.matrix(R))
      stop("Could not extract correlation matrix from gee object.",
           call. = FALSE)
    return(R)
  }

  # For all remaining classes we reconstruct from compact alpha + corstr
  if (missing(max_size) || is.null(max_size))
    stop("max_size must be supplied for ", paste(class(fit), collapse = "/"),
         " objects.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # geepack::geeglm — reconstruct from fit$geese$alpha + fit$corstr
  #
  # corstr strings used by geepack: "independence", "exchangeable", "ar1",
  # "unstructured", "userdefined", "fixed"
  # fit$geese$alpha is the compact vector from the C++ engine.
  # ---------------------------------------------------------------------------
  if (inherits(fit, "geeglm")) {
    corstr <- tolower(fit$corstr)
    alpha  <- fit$geese$alpha
    Mv     <- if (!is.null(attr(alpha, "Mv"))) attr(alpha, "Mv") else 1L
    return(.build_corr_matrix(corstr, alpha, max_size, Mv))
  }

  # ---------------------------------------------------------------------------
  # mmmgee::geem2 — reconstruct from fit$alpha + fit$corstr
  # ---------------------------------------------------------------------------
  if (inherits(fit, "geem2")) {
    corstr <- tolower(fit$corstr)
    alpha  <- fit$alpha
    Mv     <- if (!is.null(fit$Mv)) as.integer(fit$Mv) else 1L
    return(.build_corr_matrix(corstr, alpha, max_size, Mv))
  }

  stop("Unsupported class: ", paste(class(fit), collapse = ", "),
       ". Supported classes: geeglm, gee, glmgee, geem2.", call. = FALSE)
}


# =============================================================================
# Internal: build a correlation matrix from structure name + compact alpha.
# =============================================================================
.build_corr_matrix <- function(corstr, alpha, size, Mv = 1L) {

  R <- diag(size)

  if (corstr %in% c("independence", "independent")) {
    # R stays identity; alpha is unused

  } else if (corstr %in% c("exchangeable", "exch")) {
    R[lower.tri(R)] <- alpha[1L]
    R[upper.tri(R)] <- alpha[1L]

  } else if (corstr == "ar1") {
    # R[i,j] = alpha^|i-j|
    for (i in seq_len(size))
      for (j in seq_len(size))
        if (i != j) R[i, j] <- alpha[1L]^abs(i - j)

  } else if (corstr %in% c("m-dependent", "stat_m_dep", "non_stat_m_dep",
                           "ar-m")) {
    Mv <- as.integer(Mv)
    for (k in seq_len(min(Mv, size - 1L))) {
      val <- if (k <= length(alpha)) alpha[k] else 0
      for (i in seq_len(size))
        for (j in seq_len(size))
          if (abs(i - j) == k) R[i, j] <- val
    }

  } else if (corstr %in% c("unstructured", "unstr", "un")) {
    # alpha is the lower triangle in column-major order
    if (length(alpha) != size * (size - 1L) / 2L)
      stop("length(alpha) does not match choose(max_size, 2) for unstructured.",
           call. = FALSE)
    R[lower.tri(R)] <- alpha
    R[upper.tri(R)] <- t(R)[upper.tri(R)]

  } else {
    stop("Cannot reconstruct full matrix for corstr = '", corstr, "'.",
         call. = FALSE)
  }

  R
}
