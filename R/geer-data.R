#' Cerebrovascular
#'
#' The data are from a randomized crossover trial on cerebrovascular deficiency.
#'
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{period}{numeric vector indicating the period identifier.}
#'   \item{ecg}{numeric vector indicating the ecg status, coded as (0) for
#'              abnormal electrocardiogram and as (1) for normal electrocardiogram.}
#'   \item{treatment}{factor indicating the treatment group: active and placebo.}
#'   }
#'
#' @details
#' Sixty-seven subjects were enrolled in a two-period crossover trial. Subjects were
#' randomly assigned to receive either placebo followed by active treatment or active
#' treatment followed by placebo. ECG status was classified as normal or abnormal.
#'
#' @references
#' Jones, B. and Kenward, M.G. (1989) \emph{Design and Analysis of Cross-over
#' Trials.} London: Chapman and Hall/CRC Press.
#'
#' @examples
#' data("cerebrovascular")
#' str(cerebrovascular)
"cerebrovascular"


#' Cholecystectomy
#'
#' The data are from a randomized clinical trial on shoulder pain after
#' laparoscopic cholecystectomy.
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
#' Forty-one subjects were enrolled in a randomized clinical trial of shoulder pain
#' after laparoscopic cholecystectomy. Subjects were assigned to either the active
#' group (with abdominal suction) or the placebo group (without abdominal suction).
#' Shoulder pain was assessed on six occasions: morning and afternoon for three
#' days postoperatively.
#'
#' @references
#' Jorgensen J.O., Gillies R.B., Hunt D.R., Caplehorn J.R.M. and
#' Lumley T. (1995) A simple and effective way to reduce postoperative pain after
#' laparoscopic cholecystectomy. \emph{Australian and New Zealand Journal of Surgery},
#' \strong{65}, 466--469.
#'
#' @source
#' Lumley T. (1996) Generalized estimating equations for ordinal data:
#' A note on working correlation structures. \emph{Biometrics}, \strong{52}, 354--361.
#'
#' @examples
#' data("cholecystectomy")
#' str(cholecystectomy)
"cholecystectomy"


#' Depression
#'
#' The data are from a randomized clinical trial on the efficacy of oestrogen
#' patches in treating postnatal depression.
#'
#' @format A data frame with 366 rows and 5 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{visit}{numeric vector indicating the follow-up month.}
#'   \item{score}{numeric vector indicating the Edinburgh Postnatal Depression Scale
#'                (EDPS) score.}
#'   \item{treatment}{factor indicating the treatment group: placebo and oestrogen.}
#'   \item{baseline}{numeric vector indicating the EDPS score at baseline.}
#'   }
#'
#' @details
#' Sixty-one women with major depression were randomly assigned to a placebo control
#' group or an oestrogen patch group. Before treatment, depressive symptoms were
#' assessed with the EPDS. EPDS scores were then collected monthly for six months.
#' Higher scores indicate greater depression.
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
#' The data are from a randomized clinical trial on the efficacy of progabide in
#' treating partial seizures.
#'
#' @format A data frame with 236 rows and 6 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{visit}{numeric vector indicating which of the two-week interval
#'                corresponds to the reported number of epileptic seizures.}
#'   \item{seizures}{numeric vector indicating the number of epileptic seizures.}
#'   \item{treatment}{factor indicating the treatment group: placebo and progabide.}
#'   \item{lnbaseline}{numeric vector indicating the logarithm of one-quarter of
#'                     the number of epileptic seizures in the baseline 8-week
#'                     interval.}
#'   \item{lnage}{numeric vector indicating the logarithm of age in years.}
#'   }
#'
#' @details
#' Fifty-nine subjects with partial seizures were enrolled in a randomized clinical
#' trial. Subjects were assigned to receive either the anti-epileptic drug progabide
#' or placebo. Seizure counts were recorded during the 8 weeks prior to treatment.
#' After treatment, seizure counts were assessed at four 2-week clinic visits.
#'
#' @source
#' Thall, P.F. and Vail, S.C. (1990) Some covariance models for longitudinal
#' count data with overdispersion. \emph{Biometrics}, \strong{46}, 657--671.
#'
#' @references
#' Carey V.J., Wang Y.G. (2011) Working covariance model selection for generalized
#' estimating equations. \emph{Statistics in Medicine}, \strong{30}, 3117--3124.
#'
#' @examples
#' data("epilepsy")
#' str(epilepsy)
"epilepsy"


#' Leprosy
#'
#' The data are from a randomized clinical trial on the efficacy of antibiotic
#' treatments for leprosy in the Philippines.
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
#' Thirty subjects in the Philippines were enrolled in a randomized clinical trial.
#' Subjects were assigned to receive treatment A, B, or C. Before treatment, the
#' number of leprosy bacilli at six body sites was recorded. After treatment, bacilli
#' counts were recorded again, and respiratory status was assessed at baseline and at
#' four follow-up visits. The trial aimed to test whether antibiotic treatments A and B
#' reduced bacilli abundance compared with placebo (treatment C).
#'
#' @references
#' Snedecor, G.W. and Cochran W.G. (1967) \emph{Statistical Methods}.
#' Ames, Iowa: Iowa State University Press.
#'
#' @examples
#' data("leprosy")
#' str(leprosy)
"leprosy"


#' Respiratory
#'
#' The data are from a randomized clinical trial on respiratory illness.
#'
#' @format A data frame with 444 rows and 8 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{visit}{numeric vector indicating the visit identifier.}
#'   \item{status}{numeric vector indicating the respiratory status, coded as (0)
#'                 for poor and as (1) for good.}
#'   \item{treatment}{factor indicating the treatment group: active and placebo.}
#'   \item{baseline}{numeric vector indicating the respiratory status at the
#'                   baseline, coded as (0) for poor and as (1) for good.}
#'   \item{age}{numeric vector indicating the age in years recorded at baseline.}
#'   \item{gender}{factor indicating the gender: female and male.}
#'   \item{center}{factor indicating the center identifier: C1 and C2.}
#'   }
#'
#' @details
#' One hundred eleven subjects from two centers were enrolled in a randomized clinical
#' trial. Subjects were randomly assigned to treatment groups. Baseline examinations
#' were conducted before treatment. After treatment, subjects were examined at four
#' follow-up visits, and respiratory status was assessed at baseline and at each
#' follow-up.
#'
#' @references Stokes, M.E., Davis, C.S. and Koch, G.G. (1995) \emph{Categorical Data Analysis
#' using the SAS System.} Cary, NC: SAS Institute, Inc.
#'
#' @examples
#' data("respiratory")
#' str(respiratory)
"respiratory"


#' Rinse
#'
#' The data are from a randomized clinical trial on the efficacy of mouth rinses in
#' reducing dental plaque.
#'
#' @format A data frame with 210 rows and 8 columns:
#' \describe{
#'   \item{id}{numeric vector indicating the subject identifier.}
#'   \item{time}{numeric vector indicating the follow-up month.}
#'   \item{score}{numeric vector indicating the subject's score of plaque.}
#'   \item{treatment}{factor indicating the the type of rinse group: A, B and placebo.}
#'   \item{baseline}{numeric vector indicating the subject's score of plaque at
#'                   the baseline.}
#'   \item{gender}{factor indicating the gender: female and male.}
#'   \item{age}{numeric vector indicating the age in years recorded at baseline.}
#'   \item{smoke}{factor indicating the smoking status: yes and no.}
#'   }
#'
#' @details
#' One hundred nine adults aged 18–55 with pre-existing dental plaque but without
#' advanced periodontal disease were enrolled in a randomized, double-blinded
#' clinical trial. Subjects were assigned to one of two novel mouth rinses (A or B)
#' or to a control mouth rinse. Eligibility required at least 20 sound natural teeth
#' and a mean plaque index of 2.0 or greater; subjects with gross oral pathology or
#' recent antibiotic, antibacterial, or anti-inflammatory therapy were excluded.
#' Plaque was assessed at baseline, 3 months, and 6 months using the Turesky
#' modification of the Quigley-Hein index. Four subjects had missing plaque scores.
#' The trial aimed to evaluate the effectiveness of the three rinses in inhibiting
#' dental plaque.
#'
#' @references
#' Hadgu A. and Koch G. (1999) Application of generalized estimating equations
#' to a dental randomized clinical trial. \emph{Journal of Biopharmaceutical
#' Statistics},  \strong{9}, 161--178.
#'
#' @examples
#' data("rinse")
#' str(rinse)
"rinse"
