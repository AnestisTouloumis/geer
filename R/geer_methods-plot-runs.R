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
#' @param col colour of the plotting symbols. The default, \code{NULL},
#'   colours the symbols by run using \code{run_colors}. A single colour
#'   suppresses that and draws every symbol alike.
#' @param run_colors vector of colours recycled across successive runs when
#'   \code{col} is \code{NULL}. At least three are needed for adjacent runs
#'   to be distinguishable, because consecutive runs necessarily differ in
#'   sign and so are already separated vertically.
#' @param pch,cex plotting symbol and symbol expansion, passed to
#'   \code{\link[graphics]{points}}.
#' @param main,xlab,ylab character strings giving the title and the axis
#'   labels. Defaults are supplied when these are \code{NULL}.
#' @param ... further graphical parameters passed to
#'   \code{\link[graphics]{plot.default}}.
#'
#' @details
#' Each retained residual contributes one point at \eqn{+1} or \eqn{-1},
#' plotted against its position in the sequence that was tested. Residuals
#' equal to zero are excluded from that sequence, so positions index the
#' tested signs rather than the observations of the fitted model.
#'
#' Successive runs are drawn in different colours, recycled from
#' \code{run_colors}, so that a long run appears as a monochrome block and
#' the number of colour changes is the observed number of runs less one.
#' Colouring is used in preference to vertical rules at every sign change,
#' which saturate the display once the sequence runs to a few hundred
#' observations. Note that colouring by run parity alone would carry no
#' information, since consecutive runs differ in sign by definition and are
#' therefore already separated on the vertical axis; at least three colours
#' are required. Dotted grid lines mark the boundaries
#' between clusters, which makes it possible to see how many clusters have
#' residuals of a common sign; this is the pattern that drives the statistic
#' downwards under the natural ordering. Cluster boundaries are only
#' interpretable when the sequence is ordered naturally, since any other
#' ordering interleaves the clusters, and a warning is issued if they are
#' requested for such an ordering.
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
                                cex = 1L,
                                main = NULL,
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
      "cluster boundaries are only interpretable under the natural cluster/repeated ordering",
      call. = FALSE
    )
  }

  sequence_length <- length(signs)
  position <- seq_len(sequence_length)
  if (is.null(main)) {
    main <- "Residual runs"
  }
  if (is.null(xlab)) {
    xlab <- "Position in the tested sequence"
  }
  if (is.null(ylab)) {
    ylab <- "Residual sign"
  }

  graphics::plot.default(
    position,
    signs,
    type = "n",
    yaxt = "n",
    ylim = c(-1.5, 1.5),
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  graphics::axis(2L, at = c(-1, 1), labels = c("-", "+"))
  if (cluster_breaks && length(x$cluster) == sequence_length) {
    boundaries <- which(
      x$cluster[-1L] != x$cluster[-sequence_length]
    ) + 0.5
    if (length(boundaries) > 0L) {
      graphics::abline(v = boundaries, col = "grey70", lty = 3L)
    }
  }
  graphics::abline(h = 0, col = "grey40")
  ## Runs are indexed by the cumulative number of sign changes, so every
  ## element of a run receives the same colour and the colour changes exactly
  ## at the run boundaries.
  run_index <- cumsum(c(1L, as.integer(signs[-1L] != signs[-sequence_length])))
  point_colors <- if (is.null(col)) {
    if (length(run_colors) < 1L) {
      stop("'run_colors' must contain at least one colour", call. = FALSE)
    }
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
      cluster = if (length(x$cluster) == sequence_length) {
        x$cluster
      } else {
        rep(NA_real_, sequence_length)
      }
    )
  )
}
