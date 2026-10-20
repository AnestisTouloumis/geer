#' @title
#' Stepwise Model Selection by an Information Criterion
#'
#' @description
#' Performs stepwise model selection for a fitted \code{geer} object, moving at
#' each step to the single-term addition or deletion that most reduces a
#' quasi-likelihood or Gaussian pseudolikelihood information criterion.
#'
#' @inheritParams step_p
#' @param criterion character string specifying the information criterion
#'   minimized in the search: \code{"QICu"}, \code{"QIC"}, \code{"QICC"},
#'   \code{"AGPC"}, or \code{"SGPC"}. See Details. Defaults to \code{"QICu"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used in the penalty terms of \code{"QIC"} and \code{"QICC"}: one of the
#'   estimators described in \code{\link{vcov.geer}}. It is ignored for the
#'   other criteria. Defaults to \code{"robust"}, which gives the classical
#'   definitions.
#' @param trace logical indicating whether the candidate table of each step is
#'   printed during the search. Defaults to \code{FALSE}.
#'
#' @details
#' The criteria are those of \code{\link{geecriteria}} that are intended for
#' comparing mean models: \code{QIC} and \code{QICu} (Pan, 2001), \code{QICC}
#' (Hardin and Hilbe, 2013), and \code{AGPC} and \code{SGPC} (Carey and Wang,
#' 2011; Fu, Hao and Wang, 2018), as used for variable selection by
#' \pkg{glmtoolbox} (Vanegas, Rondon and Paula, 2023). The criteria that
#' \code{\link{geecriteria}} provides for selecting the working association
#' structure, such as \code{CIC}, are not accepted.
#'
#' The initial model must contain every term of \code{scope$lower}, and
#' \code{scope$upper} every term of the initial model. At each step every
#' candidate model, obtained by adding a term of \code{scope$upper} or deleting
#' a term not in \code{scope$lower}, is refitted to the rows used by
#' \code{object}, as described in \code{\link{add1.geer}}. The association and
#' dispersion parameters are re-estimated for each candidate. With
#' \code{direction = "backward"} only deletions, with \code{"forward"} only
#' additions, and with \code{"both"} additions and deletions together are
#' considered. The candidate with the smallest value of the criterion is
#' selected, the first in the order of the terms of the scope in case of a
#' tie, and the search moves to it only if its value is smaller than that of
#' the current model by more than
#' \code{sqrt(.Machine$double.eps) * max(1, abs(value))}, so that differences
#' within the convergence tolerance of the refits do not decide a move. Moves
#' respect marginality. The search stops when no candidate improves the
#' criterion or after \code{steps} steps. Because each move decreases the
#' criterion, no model is visited twice. A candidate whose criterion cannot be
#' evaluated is skipped.
#'
#' Each step refits every candidate model, but a model already fitted during
#' the search, such as the model before the previous move when
#' \code{direction = "both"}, is not refitted. With
#' \code{cov_type = "jackknife"} and \code{criterion = "QIC"} or \code{"QICC"},
#' every candidate model is in addition refitted once per cluster, which can be
#' slow.
#'
#' The quasi-likelihood in \code{QIC}, \code{QICu}, and \code{QICC} diverges
#' when a fitted mean reaches the boundary of its support, so these criteria
#' should not be used when only some of the candidate models reach the
#' boundary. Refits of such candidates warn about fitted probabilities or
#' rates near the boundary. See \code{\link{geecriteria}}.
#'
#' @return
#' The selected model, a fitted model object of class \code{"geer"} with an
#' additional component \code{anova} of class \code{c("geer_anova", "anova",
#' "data.frame")} with one row for the initial model and one for each step,
#' whose row names give the term added or removed, and columns \code{Step}
#' (\code{NA} for the initial model), \code{Df}, and the value of the
#' criterion, named after it, and an additional component \code{candidates},
#' a list with one table of class \code{c("geer_anova", "anova",
#' "data.frame")} per step. Each table has one row for the current model
#' (\code{<none>}) and one for each candidate, sorted by the criterion, and
#' columns \code{Df} and the value of the criterion. The last table is that of
#' the step at which the search stopped, unless it stopped at the step limit.
#'
#' @references
#' Carey, V.J. and Wang, Y.G. (2011) Working covariance model selection for
#' generalized estimating equations. \emph{Statistics in Medicine},
#' \bold{30}, 3117--3124.
#'
#' Fu, L., Hao, Y. and Wang, Y.G. (2018) Working correlation structure
#' selection in generalized estimating equations. \emph{Computational
#' Statistics}, \bold{33}, 983--996.
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' Pan, W. (2001) Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, \bold{57}, 120--125.
#'
#' Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
#' Equations using the new R package glmtoolbox. \emph{The R Journal},
#' \bold{15}, 105--133.
#'
#' @seealso \code{\link{step_p}}, \code{\link{geecriteria}},
#'   \code{\link{add1.geer}}, \code{\link{drop1.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' \donttest{
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fit <- geewa_binary(
#'   formula = status ~ baseline + treatment + gender + visit + age,
#'   link = "probit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   method = "pgee-jeffreys"
#' )
#' ## Backward elimination by QICu
#' step_criterion(fit)
#'
#' ## Bidirectional selection by QIC within a scope
#' step_criterion(
#'   fit,
#'   scope = list(
#'     lower = ~ baseline + treatment,
#'     upper = ~ (baseline + treatment + gender + visit + age)^2
#'   ),
#'   direction = "both",
#'   criterion = "QIC"
#' )
#' }
#'
#' @export
step_criterion <- function(object,
                           scope,
                           direction = c("backward", "forward", "both"),
                           criterion = c("QICu", "QIC", "QICC", "AGPC",
                                         "SGPC"),
                           cov_type = "robust",
                           steps = 1000,
                           trace = FALSE) {
  object <- check_geer_object(object)
  direction <- match_direction_type(direction[1L])
  criterion <- match_step_criterion(criterion[1L])
  cov_type <- match_cov_type(cov_type, default = "robust")
  steps <- check_step_count(steps)
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE", call. = FALSE)
  }
  scope_value <- if (missing(scope)) NULL else scope
  if (is.null(scope_value) && identical(direction, "forward")) {
    warning(
      paste0(
        "no 'scope' was supplied for direction = \"forward\", so the initial ",
        "model is returned unchanged"
      ),
      call. = FALSE
    )
  }
  run_step_criterion(
    object = object,
    scope = scope_value,
    direction = direction,
    criterion = criterion,
    cov_type = cov_type,
    steps = steps,
    trace = trace
  )
}
