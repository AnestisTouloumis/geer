#' Cerebrovascular
#'
#' Data from a randomized two-period crossover trial on cerebrovascular deficiency.
#'
#' @docType data
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{period}{integer period identifier.}
#'   \item{ecg}{integer indicator for ECG status, 0 = abnormal, 1 = normal.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#' }
#'
#' @details
#' Sixty-seven subjects were enrolled in a two-period crossover trial. Subjects were
#' randomly assigned to receive either placebo followed by active treatment or active
#' treatment followed by placebo. ECG status was classified as normal or abnormal.
#'
#' @references
#' Jones, B. and Kenward, M.G. (1989) \emph{Design and Analysis of Cross-over Trials.}
#' London: Chapman and Hall/CRC Press.
#'
#' @examples
#' data("cerebrovascular")
#' str(cerebrovascular)
"cerebrovascular"


#' Cholecystectomy
#'
#' Data from a randomized clinical trial on shoulder pain after laparoscopic cholecystectomy.
#'
#' @docType data
#' @format A data frame with 246 rows and 6 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{time}{integer time identifier.}
#'   \item{pain}{integer indicator for shoulder pain status, 0 = high pain, 1 = low pain.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{age}{numeric age in years.}
#' }
#'
#' @details
#' Forty-one subjects were enrolled in a randomized clinical trial of shoulder pain
#' after laparoscopic cholecystectomy. Subjects were assigned to either the active
#' group (with abdominal suction) or the placebo group (without abdominal suction).
#' Shoulder pain was assessed on six occasions: morning and afternoon for three
#' days postoperatively.
#'
#' @references
#' Jorgensen J.O., Gillies R.B., Hunt D.R., Caplehorn J.R.M. and Lumley T. (1995)
#' A simple and effective way to reduce postoperative pain after laparoscopic cholecystectomy.
#' \emph{Australian and New Zealand Journal of Surgery}, \strong{65}, 466--469.
#'
#' @source
#' Lumley T. (1996) Generalized estimating equations for ordinal data: A note on working
#' correlation structures. \emph{Biometrics}, \strong{52}, 354--361.
#'
#' @examples
#' data("cholecystectomy")
#' str(cholecystectomy)
"cholecystectomy"


#' Depression
#'
#' Data from a randomized clinical trial on the efficacy of oestrogen patches in treating
#' postnatal depression.
#'
#' @docType data
#' @format A data frame with 366 rows and 5 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer follow-up month identifier.}
#'   \item{score}{numeric Edinburgh Postnatal Depression Scale (EPDS) score. Higher scores indicate greater depression.}
#'   \item{treatment}{factor with levels \code{placebo} and \code{oestrogen}.}
#'   \item{baseline}{numeric EPDS score at baseline.}
#' }
#'
#' @details
#' Sixty-one women with major depression were randomly assigned to a placebo control group
#' or an oestrogen patch group. Before treatment, depressive symptoms were assessed with
#' the EPDS. EPDS scores were then collected monthly for six months.
#'
#' @source
#' \url{https://stats.oarc.ucla.edu/spss/library/spss-librarypanel-data-analysis-using-gee/}
#'
#' @references
#' Gregoire A.J.P., Kumar R., Everitt B., Henderson A.F. and Studd J.W.W. (1996)
#' Transdermal oestrogen for treatment of severe postnatal depression.
#' \emph{The Lancet}, \strong{347}, 930--933.
#'
#' @examples
#' data("depression")
#' str(depression)
"depression"


#' Epilepsy
#'
#' Data from a randomized clinical trial on the efficacy of progabide in treating partial seizures.
#'
#' @docType data
#' @format A data frame with 236 rows and 6 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer visit identifier (two-week interval).}
#'   \item{seizures}{integer count of epileptic seizures.}
#'   \item{treatment}{factor with levels \code{placebo} and \code{progabide}.}
#'   \item{lnbaseline}{numeric logarithm of one-quarter of the number of seizures in the baseline 8-week interval.}
#'   \item{lnage}{numeric logarithm of age in years.}
#' }
#'
#' @details
#' Fifty-nine subjects with partial seizures were enrolled in a randomized clinical trial.
#' Subjects were assigned to receive either the anti-epileptic drug progabide or placebo.
#' Seizure counts were recorded during the 8 weeks prior to treatment. After treatment,
#' seizure counts were assessed at four 2-week clinic visits.
#'
#' @source
#' Thall, P.F. and Vail, S.C. (1990) Some covariance models for longitudinal count data with
#' overdispersion. \emph{Biometrics}, \strong{46}, 657--671.
#'
#' @references
#' Carey V.J. and Wang Y.G. (2011) Working covariance model selection for generalized estimating
#' equations. \emph{Statistics in Medicine}, \strong{30}, 3117--3124.
#'
#' @examples
#' data("epilepsy")
#' str(epilepsy)
"epilepsy"


#' Leprosy
#'
#' Data from a randomized clinical trial on the efficacy of antibiotic treatments for leprosy
#' in the Philippines.
#'
#' @docType data
#' @format A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{period}{factor period identifier with levels \code{pre} and \code{post}.}
#'   \item{bacilli}{integer count of leprosy bacilli at six sites of the body.}
#'   \item{treatment}{factor with levels \code{A}, \code{B} and \code{C} (C = placebo).}
#' }
#'
#' @details
#' Thirty subjects in the Philippines were enrolled in a randomized clinical trial. Subjects
#' were assigned to receive treatment A, B, or C. Before treatment, the number of leprosy
#' bacilli at six body sites was recorded. After treatment, bacilli counts were recorded again.
#' The trial aimed to test whether antibiotic treatments A and B reduced bacilli abundance
#' compared with placebo (treatment C).
#'
#' @references
#' Snedecor, G.W. and Cochran W.G. (1967) \emph{Statistical Methods}. Ames, Iowa: Iowa State University Press.
#'
#' @examples
#' data("leprosy")
#' str(leprosy)
"leprosy"


#' Respiratory
#'
#' Data from a randomized clinical trial on respiratory illness.
#'
#' @docType data
#' @format A data frame with 444 rows and 8 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{visit}{integer visit identifier.}
#'   \item{status}{integer indicator for respiratory status, 0 = poor, 1 = good.}
#'   \item{treatment}{factor with levels \code{active} and \code{placebo}.}
#'   \item{baseline}{integer indicator for baseline respiratory status, 0 = poor, 1 = good.}
#'   \item{age}{numeric age in years recorded at baseline.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{center}{factor with levels \code{C1} and \code{C2}.}
#' }
#'
#' @details
#' One hundred eleven subjects from two centers were enrolled in a randomized clinical trial.
#' Subjects were randomly assigned to treatment groups. Baseline examinations were conducted
#' before treatment. After treatment, subjects were examined at four follow-up visits.
#'
#' @references
#' Stokes, M.E., Davis, C.S. and Koch, G.G. (1995) \emph{Categorical Data Analysis using the SAS System.}
#' Cary, NC: SAS Institute, Inc.
#'
#' @examples
#' data("respiratory")
#' fit <- geewa_binary(status ~ treatment + baseline + age + gender,
#'                     id = id, data = respiratory,
#'                     link = "logit",
#'                     orstr = "exchangeable")
#' summary(fit)
"respiratory"


#' Rinse
#'
#' Data from a randomized clinical trial on the efficacy of mouth rinses in reducing dental plaque.
#'
#' @docType data
#' @format A data frame with 218 rows and 8 columns:
#' \describe{
#'   \item{id}{integer subject identifier.}
#'   \item{time}{integer follow-up month identifier.}
#'   \item{score}{numeric plaque score.}
#'   \item{treatment}{factor indicating the type of rinse group with levels \code{A}, \code{B} and \code{placebo}.}
#'   \item{baseline}{numeric plaque score at baseline.}
#'   \item{gender}{factor with levels \code{female} and \code{male}.}
#'   \item{age}{numeric age in years recorded at baseline.}
#'   \item{smoke}{factor with levels \code{yes} and \code{no}.}
#' }
#'
#' @details
#' One hundred nine adults aged 18-55 with pre-existing dental plaque but without advanced
#' periodontal disease were enrolled in a randomized, double-blinded clinical trial. Subjects
#' were assigned to one of two novel mouth rinses (A or B) or to a control mouth rinse.
#' Eligibility required at least 20 sound natural teeth and a mean plaque index of 2.0 or greater.
#' Plaque was assessed at baseline, 3 months, and 6 months using the Turesky modification of the
#' Quigley-Hein index. Four subjects had missing plaque scores. The trial aimed to evaluate the
#' effectiveness of the three rinses in inhibiting dental plaque.
#'
#' @references
#' Hadgu A. and Koch G. (1999) Application of generalized estimating equations to a dental randomized clinical trial.
#' \emph{Journal of Biopharmaceutical Statistics}, \strong{9}, 161--178.
#'
#' @examples
#' data("rinse")
#' str(rinse)
"rinse"
