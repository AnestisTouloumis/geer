#' @title
#' Plot the Residual Sign Sequence of a Runs Test
#'
#' @method plot geer_runs_test
#'
#' @description
#' Plots the sequence of residual signs tested by \code{\link{runs_test}}, as
#' in Hardin and Hilbe (2013, Section 4.2.1), with successive runs in
#' different colors and, under the natural ordering, the cluster boundaries.
#'
#' @param x an object of class \code{"geer_runs_test"}, as returned by
#'   \code{\link{runs_test}}.
#' @param cluster_breaks logical indicating whether dotted lines mark the
#'   boundaries between clusters. Defaults to \code{NULL}, which draws them
#'   only under the natural ordering.
#' @param col color of the symbols. Defaults to \code{NULL}, which colors the
#'   symbols by run using \code{run_colors}.
#' @param run_colors vector of at least three colors, recycled across
#'   successive runs. Defaults to
#'   \code{grDevices::palette.colors(4L, "Okabe-Ito")}.
#' @param pch plotting symbol, passed to \code{\link[graphics]{points}}.
#'   Defaults to \code{1L}.
#' @param cex symbol size, passed to \code{\link[graphics]{points}}. Defaults to
#'   \code{NULL}, which shrinks the symbols for long sequences.
#' @param main,sub,xlab,ylab character strings giving the title, subtitle, and
#'   axis labels. Defaults to \code{NULL}, which gives informative labels; use
#'   \code{""} to suppress one.
#' @param ... further named graphical parameters passed to
#'   \code{\link[graphics]{plot.default}}, such as \code{ylim}. \code{type}
#'   cannot be changed.
#'
#' @details
#' Each tested residual is plotted at \eqn{+1} or \eqn{-1} against its
#' position in the tested sequence; zero residuals are excluded and their
#' number is given in the default subtitle. The default title gives the
#' ordering, and the default subtitle the number of runs, its null
#' expectation, and the p-value. At least three colors are needed, because
#' with two the color would simply repeat the sign. Cluster boundaries are
#' meaningful only under the natural ordering, and requesting them otherwise
#' gives a warning.
#'
#' @return
#' Invisibly, a data frame with one row per tested residual and columns
#' \code{position}, \code{sign}, \code{run}, and \code{cluster}.
#'
#' @references
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' @seealso \code{\link{runs_test}}, \code{\link{residuals.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   repeated = visit,
#'   corstr = "exchangeable"
#' )
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
        "'cluster_breaks' must be NULL, TRUE or FALSE",
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
        "'run_colors' must contain at least three colors; with two, the ",
        "colors only repeat the signs"
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
        "'type' cannot be changed: the symbols are drawn separately to color ",
        "them by run"
      ),
      call. = FALSE
    )
  }
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
