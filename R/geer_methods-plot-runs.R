#' @title
#' Plot the Residual Sign Sequence of a Runs Test
#'
#' @method plot geer_runs_test
#'
#' @description
#' Draws the graphical companion of the Wald-Wolfowitz runs test described by
#' Hardin and Hilbe (2013, Section 4.2.1): the sequence of residual signs that
#' was tested, with the boundaries between runs and, where meaningful, the
#' boundaries between clusters.
#'
#' @param x an object of class \code{"geer_runs_test"}, as returned by
#'   \code{\link{runs_test}}.
#' @param cluster_breaks logical value indicating whether vertical grid lines
#'   should mark the boundaries between clusters. The default, \code{NULL},
#'   draws them when the test used the natural cluster/repeated ordering and
#'   omits them otherwise.
#' @param col color of the plotting symbols. The default, \code{NULL},
#'   colors the symbols by run using \code{run_colors}. A single color
#'   suppresses that and draws every symbol alike; a vector is recycled across
#'   the tested residuals by \code{\link[graphics]{points}}.
#' @param run_colors vector of at least three colors recycled across
#'   successive runs when \code{col} is \code{NULL}. Fewer than three is an
#'   error; see Details.
#' @param pch plotting symbol, passed to \code{\link[graphics]{points}}.
#' @param cex symbol expansion, passed to \code{\link[graphics]{points}}. The
#'   default, \code{NULL}, shrinks the symbols as the sequence lengthens so
#'   that long sequences remain readable.
#' @param main,sub,xlab,ylab character strings giving the title, subtitle and
#'   axis labels. Defaults are supplied when these are \code{NULL}; use
#'   \code{""} to suppress one.
#' @param ... further graphical parameters passed to
#'   \code{\link[graphics]{plot.default}}. These take precedence over the
#'   defaults set here, so \code{ylim} and \code{yaxt} may be overridden;
#'   \code{type} may not, because the symbols are drawn separately. All must
#'   be named.
#'
#' @details
#' Each retained residual contributes one point at \eqn{+1} or \eqn{-1},
#' plotted against its position in the sequence that was tested. Residuals
#' equal to zero are excluded from that sequence, so positions index the
#' tested signs rather than the observations of the fitted model. The number
#' of excluded residuals is reported in the default subtitle when it is not
#' zero.
#'
#' The default title records the ordering that was used, because the ordering
#' determines which alternative the test has power against, and two sequences
#' drawn from the same fit under different orderings are otherwise
#' indistinguishable. The default subtitle reports the observed number of
#' runs, its null expectation and the p-value, so that the figure can be read
#' without the printed test beside it.
#'
#' Successive runs are drawn in different colors, recycled from
#' \code{run_colors}, so that a long run appears as a monochrome block and
#' the number of color changes is the observed number of runs less one.
#' Coloring is used in preference to vertical rules at every sign change,
#' which saturate the display once the sequence runs to a few hundred
#' observations. At least three colors are required: with two, the color
#' alternates exactly with the sign, so it reproduces the vertical axis and
#' carries no information of its own.
#'
#' Dotted grid lines mark the boundaries between clusters, which makes it
#' possible to see how many clusters have residuals of a common sign; this is
#' the pattern that drives the statistic downwards under the natural ordering.
#' Cluster boundaries are only interpretable when the sequence is ordered
#' naturally, since any other ordering interleaves the clusters, and a warning
#' is issued if they are requested for such an ordering.
#'
#' @return
#' Invisibly, a data frame with one row per tested residual and columns
#' \code{position}, \code{sign}, \code{run} and \code{cluster}, where
#' \code{run} numbers the runs consecutively from one.
#'
#' @references
#' Chang, Y.-C. (2000) Residuals analysis of the generalized linear models for
#' longitudinal data. \emph{Statistics in Medicine}, \bold{19}, 1277--1293.
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' @seealso \code{\link{runs_test}}, \code{\link{residuals.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   seizures ~ treatment + lnbaseline + lnage,
#'   data = epilepsy,
#'   id = id,
#'   repeated = visit,
#'   family = poisson(link = "log"),
#'   corstr = "exchangeable"
#' )
#'
#' plot(runs_test(fit))
#' plot(runs_test(fit, order_by = "fitted"), ylim = c(-2, 2))
#'
#' @export
plot.geer_runs_test <- function(x,
                                cluster_breaks = NULL,
                                col = NULL,
                                run_colors = grDevices::palette.colors(
                                  4L,
                                  "Okabe-Ito"
                                ),
                                pch = 1L,
                                cex = NULL,
                                main = NULL,
                                sub = NULL,
                                xlab = NULL,
                                ylab = NULL,
                                ...) {
  if (!inherits(x, "geer_runs_test")) {
    stop("'x' must be a 'geer_runs_test' object", call. = FALSE)
  }
  signs <- x$signs
  if (!is.numeric(signs) || length(signs) < 2L) {
    stop(
      "'x' must contain the tested sign sequence: refit with runs_test()",
      call. = FALSE
    )
  }
  sequence_length <- length(signs)
  cluster <- x$cluster
  if (length(cluster) != sequence_length) {
    stop(
      paste0(
        "'x' is malformed: 'cluster' must have one value per tested ",
        "residual"
      ),
      call. = FALSE
    )
  }
  if (!is.null(cluster_breaks)) {
    if (length(cluster_breaks) != 1L || !is.logical(cluster_breaks) ||
        is.na(cluster_breaks)) {
      stop(
        "'cluster_breaks' must be NULL or a single non-missing logical value",
        call. = FALSE
      )
    }
  }
  natural <- isTRUE(x$natural_order)
  if (is.null(cluster_breaks)) {
    cluster_breaks <- natural
  } else if (cluster_breaks && !natural) {
    warning(
      paste0(
        "cluster boundaries are only interpretable under the natural ",
        "cluster/repeated ordering"
      ),
      call. = FALSE
    )
  }
  if (is.null(col) && length(run_colors) < 3L) {
    stop(
      paste0(
        "'run_colors' must contain at least three colors: with two the ",
        "color alternates with the residual sign and duplicates the ",
        "vertical axis"
      ),
      call. = FALSE
    )
  }

  position <- seq_len(sequence_length)
  if (is.null(main)) {
    main <- sprintf("Residual runs ordered by %s", x$order_by)
  }
  if (is.null(sub)) {
    sub <- sprintf(
      "T = %d, E(T) = %.4g, p = %.3g",
      x$runs,
      x$expected_runs,
      x$p.value
    )
    if (isTRUE(x$zero > 0L)) {
      sub <- sprintf("%s; %d zero residuals omitted", sub, x$zero)
    }
  }
  if (is.null(xlab)) {
    xlab <- "Position in the tested sequence"
  }
  if (is.null(ylab)) {
    ylab <- "Residual sign"
  }
  if (is.null(cex)) {
    ## Open symbols merge into a solid band once the sequence is long, which
    ## is the regime the run coloring is meant to survive.
    cex <- max(0.25, min(1, 200 / sequence_length))
  }

  dots <- list(...)
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || !all(nzchar(dot_names)))) {
    stop("arguments in '...' must all be named", call. = FALSE)
  }
  if ("type" %in% dot_names) {
    stop(
      paste0(
        "'type' must not be supplied: the symbols are drawn separately so ",
        "that they can be colored by run"
      ),
      call. = FALSE
    )
  }
  ## The defaults are assembled first and then overwritten by '...', so a
  ## user-supplied 'ylim' or 'yaxt' replaces the default instead of reaching
  ## plot.default twice.
  plot_arguments <- list(
    x = position,
    y = signs,
    type = "n",
    yaxt = "n",
    ylim = c(-1.5, 1.5),
    main = main,
    sub = sub,
    xlab = xlab,
    ylab = ylab
  )
  plot_arguments[dot_names] <- dots
  do.call(graphics::plot.default, plot_arguments)
  if (identical(plot_arguments$yaxt, "n")) {
    graphics::axis(2L, at = c(-1, 1), labels = c("-", "+"))
  }
  if (cluster_breaks) {
    boundaries <- which(
      cluster[-1L] != cluster[-sequence_length]
    ) + 0.5
    if (length(boundaries) > 0L) {
      graphics::abline(v = boundaries, col = "grey70", lty = 3L)
    }
  }
  graphics::abline(h = 0, col = "grey40")
  ## Runs are indexed by the cumulative number of sign changes, so every
  ## element of a run receives the same color and the color changes exactly
  ## at the run boundaries.
  run_index <- cumsum(c(1L, as.integer(signs[-1L] != signs[-sequence_length])))
  point_colors <- if (is.null(col)) {
    run_colors[(run_index - 1L) %% length(run_colors) + 1L]
  } else {
    col
  }
  graphics::points(position, signs, pch = pch, cex = cex, col = point_colors)

  invisible(
    data.frame(
      position = position,
      sign = signs,
      run = run_index,
      cluster = cluster
    )
  )
}
