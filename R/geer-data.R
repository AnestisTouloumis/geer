#' Cerebrovascular
#'
#' The data are from a crossover trial on cerebrovascular deficiency.
#'
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{period}{numeric vector indicating the period identifier.}
#'   \item{ecg}{numeric vector indicating the ecg status, coded as (0) for
#'   abnormal electrocardiogram and as (1) for normal electrocardiogram.}
#'   \item{treatment}{factor indicating the treatment group: active and placebo.}
#'   }
#'
#' @details
#' A total of 67 subjects were enrolled in a two-period crossover trial. Of those,
#' 34 subjects received placebo followed by active treatment and 33 received
#' active treatment followed by placebo. The ecg status was determined on whether
#' the electrocardiogram was abnormal or normal.
#'
#' @references Jones, B. and Kenward, M.G. (1989) \emph{Design and Analysis of Cross-over
#' Trials.} London: Chapman and Hall/CRC Press.
#'
#' @examples
#' data("cerebrovascular")
#' str(cerebrovascular)
"cerebrovascular"

#' Cholecystectomy
#'
#' The data are from a clinical trial on shoulder pain after laparoscopic
#' cholecystectomy.
#'
#' @format A data frame with 246 rows and 6 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{time}{numeric vector indicating the time identifier.}
#'  \item{pain}{numeric vector indicating the self-reported shoulder pain status,
#'              coded as (0) for high pain and as (1) for low pain.}
#'  \item{treatment}{factor indicating the treatment group: active and placebo.}
#'   \item{gender}{factor indicating the gender: female and male.}
#'   \item{age}{numeric vector indicating the age in years.}
#'   }
#'
#' @details
#' A total of 41 subjects were enrolled in a randomized clinical trial of shoulder
#' pain after laparoscopic cholecystectomy. Subjects were allocated to either
#' the active group (with abdominal suction) or to the placebo group (without
#' abdominal suction). After the operation, subjects were asked to rate their
#' shoulder pain in six occasions (morning and afternoon for three days after
#' the operation).
#'
#' @references Jorgensen J.O., Gillies R.B., Hunt D.R., Caplehorn J.R.M. and
#' Lumley T. (1995) A simple and effective way to reduce postoperative pain after
#' laparoscopic cholecystectomy. \emph{Australian and New Zealand Journal of Surgery},
#' \strong{65}, 466-–469.
#'
#' @source Lumley T. (1996) Generalized estimating equations for ordinal data:
#' A note on working correlation structures. \emph{Biometrics}, \strong{52}, 354-–361.
#'
#' @examples
#' data("cholecystectomy")
#' str(cholecystectomy)
"cholecystectomy"


#' Epilepsy
#'
#' The data are from a clinical trial on epilepsy.
#'
#' @format A data frame with 236 rows and 6 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{visit}{numeric vector indicating which of the two-week interval
#'                corresponds to the reported number of epileptic seizures.}
#'   \item{seizures}{numeric vector indicating the number of epileptic seizures.}
#'   \item{treatment}{factor indicating the treatment group: placebo and progabide.}
#'   \item{base}{numeric vector indicating the number of epileptic seizures in the
#'               baseline 8-week interval.}
#'   \item{age}{numeric vector indicating the age in years.}
#'   }
#'
#' @details
#' A total of 59 subjects with partial seizures were enrolled in a randomized
#' clinical trial. Subjects were allocated to either the anti-epileptic drug
#' progabide or to placebo. Prior to receiving treatment, the number of epileptic
#' seizures during the preceding 8-week interval were recorded. After receiving
#' treatment, the number of epileptic seizures were recorded at each of four
#' 2-week intervals clinic visits.
#'
#' @source
#' Thall, P.F. and Vail, S.C. (1990) Some covariance models for longitudinal
#' count data with overdispersion. \emph{Biometrics}, \strong{46}, 657--671.
#'
#' @references
#' Carey V.J., Wang Y.G. (2011) Working covariance model selection for generalized
#' estimating equations. \emph{Statistics in Medicine}, \strong{30}, 3117-–3124.
#'
#' @examples
#' data("epilepsy")
#' str(epilepsy)
"epilepsy"

#' Leprosy
#'
#' The data are from a clinical trial on leprosy.
#'
#' @format A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{period}{factor indicating the period identifier: pre and post.}
#'   \item{bacilli}{numeric vector indicating the number of leprosy bacilli
#'                  at six sites of the body.}
#'   \item{treatment}{factor indicating the treatment group: A, B and C.}
#'   }
#'
#' @details
#' A total of 30 subjects were enrolled in a randomized clinical trial in
#' Philippines. Prior to receiving treatment, the number of leprosy bacilli at
#' six sites of the body were recorded. After receiving treatment (A, B or C),
#' the number of leprosy bacilli were recorded again. The respiratory status
#' was determined at the baseline and at each of the four follow-up visits. One
#' of the goals of this study is whether treatment with antibiotics
#' (treatments A and B) reduces the abundance of leprosy bacilli at the six sites
#' of the body when compared to placebo (treatment C).
#'
#' @references Snedecor, G.W. and Cochran W.G. (1967) \emph{Statistical Methods}.
#' Ames, Iowa: Iowa State University Press.
#'
#' @examples
#' data("leprosy")
#' str(leprosy)
"leprosy"


#' Respiratory
#'
#' The data are from a clinical trial on respiratory illness.
#'
#' @format A data frame with 444 rows and 8 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{visit}{numeric vector indicating the visit identifier.}
#'   \item{y}{numeric vector indicating the respiratory status, coded as (0) for
#'            poor and as (1) for good.}
#'   \item{treatment}{factor indicating the treatment group: active and placebo.}
#'   \item{baseline}{numeric vector indicating the respiratory status at the
#'                   baseline, coded as (0) for poor and as (1) for good.}
#'   \item{age}{numeric vector indicating the age in years recorded at baseline.}
#'   \item{gender}{factor indicating the gender: female and male.}
#'   \item{center}{factor indicating the center identifier: C1 and C2.}
#'   }
#'
#' @details
#' A total of 111 subjects from two different centers were enrolled in a randomized
#' clinical trial. Prior to receiving treatment, subjects were examined at baseline.
#' After receiving treatment, subjects were examined at four visits during. The
#' respiratory status was determined at the baseline and at each of the four
#' follow-up visits.
#'
#' @references Stokes, M.E., Davis, C.S. and Koch, G.G. (1995) \emph{Categorical Data Analysis
#' using the SAS System.} Cary, NC: SAS Institute, Inc.
#'
#' @examples
#' data("respiratory")
#' str(respiratory)
"respiratory"




