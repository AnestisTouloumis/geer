#' Cerebrovascular Deficiency Trial
#'
#' Data from a randomized two-period crossover trial on cerebrovascular
#' deficiency.
#'
#' @docType data
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{period}{integer period identifier.}
#'   \item{ecg}{integer indicator of ECG status, 0 = abnormal and 1 = normal.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#' }
#'
#' @details
#' Sixty-seven subjects were enrolled in a two-period crossover trial, yielding
#' 134 observations (two per subject). Subjects were randomly assigned to
#' receive either placebo followed by active treatment, or active treatment
#' followed by placebo. ECG status was classified as normal or abnormal at each
#' period.
#'
#' @references
#' Jones, B. and Kenward, M.G. (1989) \emph{Design and Analysis of Cross-over
#' Trials}. London: Chapman and Hall/CRC Press.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' str(cerebrovascular)
"cerebrovascular"


#' Shoulder Pain After Laparoscopic Cholecystectomy Trial
#'
#' Data from a randomized clinical trial on shoulder pain after laparoscopic
#' cholecystectomy.
#'
#' @docType data
#' @format A data frame with 246 rows and 6 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{time}{integer time identifier.}
#'   \item{pain}{integer indicator of shoulder pain status, 0 = high pain and
#'     1 = low pain.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{age}{numeric age in years.}
#' }
#'
#' @details
#' Forty-one subjects were enrolled in a randomized clinical trial of shoulder
#' pain after laparoscopic cholecystectomy, yielding 246 observations
#' (six per subject). Subjects were assigned to either the active group (with
#' abdominal suction) or the placebo group (without abdominal suction).
#' Shoulder pain was assessed on six occasions: morning and afternoon on each
#' of the first three postoperative days.
#'
#' @source
#' Lumley, T. (1996) Generalized estimating equations for ordinal data: a note
#' on working correlation structures. \emph{Biometrics}, \bold{52}, 354--361.
#'
#' @references
#' Jorgensen, J.O., Gillies, R.B., Hunt, D.R., Caplehorn, J.R.M. and Lumley, T.
#' (1995) A simple and effective way to reduce postoperative pain after
#' laparoscopic cholecystectomy. \emph{Australian and New Zealand Journal of
#' Surgery}, \bold{65}, 466--469.
#'
#' @examples
#' data("cholecystectomy", package = "geer")
#' str(cholecystectomy)
"cholecystectomy"


#' Postnatal Depression Oestrogen Patch Trial
#'
#' Data from a randomized clinical trial on the efficacy of oestrogen patches
#' in treating postnatal depression.
#'
#' @docType data
#' @format A data frame with 366 rows and 5 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer follow-up month identifier.}
#'   \item{score}{numeric Edinburgh Postnatal Depression Scale (EPDS) score.
#'     Higher scores indicate greater depression.}
#'   \item{treatment}{factor with levels \code{placebo} and \code{oestrogen}.}
#'   \item{baseline}{numeric EPDS score at baseline.}
#' }
#'
#' @details
#' Sixty-one women with major depression were randomly assigned to either a
#' placebo control group or an oestrogen patch group, yielding 366
#' observations (six monthly assessments per subject). Before treatment,
#' depressive symptoms were assessed using the Edinburgh Postnatal Depression
#' Scale (EPDS). EPDS scores were then collected monthly for six months.
#'
#' @source
#' \url{https://stats.oarc.ucla.edu/spss/library/spss-librarypanel-data-analysis-using-gee/}
#'
#' @references
#' Gregoire, A.J.P., Kumar, R., Everitt, B., Henderson, A.F. and Studd, J.W.W.
#' (1996) Transdermal oestrogen for treatment of severe postnatal depression.
#' \emph{The Lancet}, \bold{347}, 930--933.
#'
#' @examples
#' data("depression", package = "geer")
#' str(depression)
"depression"


#' Progabide Epilepsy Trial
#'
#' Data from a randomized clinical trial on the efficacy of progabide in
#' treating partial seizures.
#'
#' @docType data
#' @format A data frame with 236 rows and 6 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer visit identifier corresponding to a two-week interval.}
#'   \item{seizures}{integer count of epileptic seizures.}
#'   \item{treatment}{factor with levels \code{placebo} and \code{progabide}.}
#'   \item{lnbaseline}{numeric logarithm of one-quarter of the number of
#'     seizures in the baseline 8-week interval.}
#'   \item{lnage}{numeric logarithm of age in years.}
#' }
#'
#' @details
#' Fifty-nine subjects with partial seizures were enrolled in a randomized
#' clinical trial, yielding 236 observations (four visits per subject).
#' Subjects were assigned to receive either progabide or placebo. Seizure
#' counts were recorded during the 8 weeks before treatment (baseline). After
#' treatment, seizure counts were assessed at four 2-week clinic visits.
#'
#' @source
#' Thall, P.F. and Vail, S.C. (1990) Some covariance models for longitudinal
#' count data with overdispersion. \emph{Biometrics}, \bold{46}, 657--671.
#'
#' @references
#' Carey, V.J. and Wang, Y.G. (2011) Working covariance model selection for
#' generalized estimating equations. \emph{Statistics in Medicine}, \bold{30},
#' 3117--3124.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' str(epilepsy)
"epilepsy"


#' Antibiotic Treatment for Leprosy Trial
#'
#' Data from a randomized clinical trial on the efficacy of antibiotic
#' treatments for leprosy in the Philippines.
#'
#' @docType data
#' @format A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{period}{factor period identifier with levels \code{pre} and
#'     \code{post}.}
#'   \item{bacilli}{integer count of leprosy bacilli at six body sites.}
#'   \item{treatment}{factor with levels \code{A}, \code{B}, and \code{C},
#'     where \code{C} denotes placebo.}
#' }
#'
#' @details
#' Thirty subjects in the Philippines were enrolled in a randomized clinical
#' trial, yielding 60 observations (two periods per subject). Subjects were
#' assigned to receive treatment A, B, or C (placebo). Before treatment, the
#' number of leprosy bacilli at six body sites was recorded. After treatment,
#' bacilli counts were recorded again. The trial aimed to assess whether
#' treatments A and B reduced bacilli abundance compared with placebo.
#'
#' @references
#' Snedecor, G.W. and Cochran, W.G. (1967) \emph{Statistical Methods}. Ames,
#' Iowa: Iowa State University Press.
#'
#' @examples
#' data("leprosy", package = "geer")
#' str(leprosy)
"leprosy"


#' Respiratory Illness Clinical Trial
#'
#' Data from a randomized clinical trial on respiratory illness.
#'
#' @docType data
#' @format A data frame with 444 rows and 8 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer follow-up visit identifier.}
#'   \item{status}{integer indicator of respiratory status, 0 = poor and
#'     1 = good.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#'   \item{baseline}{integer indicator of baseline respiratory status,
#'     0 = poor and 1 = good.}
#'   \item{age}{numeric age in years recorded at baseline.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{center}{factor with levels \code{C1} and \code{C2}.}
#' }
#'
#' @details
#' One hundred eleven subjects from two clinical centers were enrolled in a
#' randomized clinical trial, yielding 444 observations (four follow-up
#' visits per subject). Subjects were randomly assigned to treatment groups.
#' Baseline examinations were conducted before treatment began. After
#' treatment, subjects were examined at four scheduled follow-up visits.
#'
#' @references
#' Stokes, M.E., Davis, C.S. and Koch, G.G. (1995) \emph{Categorical Data
#' Analysis using the SAS System}. Cary, NC: SAS Institute, Inc.
#'
#' @examples
#' data("respiratory", package = "geer")
#' str(respiratory)
"respiratory"


#' Dental Plaque Mouth Rinse Trial
#'
#' Data from a randomized clinical trial on the efficacy of mouth rinses in
#' reducing dental plaque.
#'
#' @docType data
#' @format A data frame with 218 rows and 8 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{time}{integer follow-up month identifier.}
#'   \item{score}{numeric plaque score.}
#'   \item{treatment}{factor indicating rinse group with levels \code{A},
#'     \code{B}, and \code{placebo}.}
#'   \item{baseline}{numeric plaque score at baseline.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{age}{numeric age in years recorded at baseline.}
#'   \item{smoke}{factor with levels \code{yes} and \code{no}.}
#' }
#'
#' @details
#' One hundred nine adults aged 18--55 with pre-existing dental plaque but
#' without advanced periodontal disease were enrolled in a randomized,
#' double-blind clinical trial, yielding 218 observations (baseline plus
#' up to two follow-up assessments per subject). Subjects were assigned to one
#' of two novel mouth rinses (A or B) or to a control mouth rinse. Eligibility
#' required at least 20 sound natural teeth and a mean plaque index of 2.0 or
#' greater. Plaque was assessed at baseline, 3 months, and 6 months using the
#' Turesky modification of the Quigley-Hein index. Four subjects had missing
#' plaque scores. The trial aimed to evaluate the effectiveness of the three
#' rinses in inhibiting dental plaque.
#'
#' @references
#' Hadgu, A. and Koch, G. (1999) Application of generalized estimating
#' equations to a dental randomized clinical trial. \emph{Journal of
#' Biopharmaceutical Statistics}, \bold{9}, 161--178.
#'
#' @examples
#' data("rinse", package = "geer")
#' str(rinse)
"rinse"
