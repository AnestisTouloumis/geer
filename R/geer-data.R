#' Two-period crossover trial study on cerebrovascular deficiency
#'
#' The dataset consists of safety data from a crossover trial on the disease
#' cerebrovascular deficiency.
#'
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{a factor with 67 levels}
#'   \item{period}{a factor with levels 1 and 2}
#'   \item{ecg}{a factor (ECG response) with levels normal and abnormal}
#'   \item{treatment}{a factor with levels placebo and active}
#'   }
#' @details
#' The response variable is not a trial endpoint but rather a potential side
#' effect. In this two-period crossover trial, comparing the effects of active
#' drug to placebo, 67 patients were randomly allocated to the two treatment
#' sequences, with 34 patients receiving placebo followed by active treatment,
#' and 33 patients receiving active treatment followed by placebo. The response
#' variable is binary, indicating whether an electrocardiogram (ECG) was
#' abnormal (1) or normal (0). Each patient has a bivariate binary response
#' vector.
#'
#' @references Jones B., Kenward M.G. (1989) Design and Analysis of Cross-over
#' Trials. London: Chapman and Hall/CRC Press.
"cerebrovascular"


#' Seizure counts for 59 epileptics
#'
#' The data are from a placebo-controlled clinical trial of 59 epileptics.
#' Patients with partial seizures were enrolled in a randomized clinical trial
#' of the anti-epileptic drug, progabide.
#'
#' @format A data frame with 295 rows and 5 columns:
#' \describe{
#'   \item{id}{a factor with 59 levels}
#'   \item{treatment}{a factor with levels placebo and progabide}
#'   \item{age}{a numeric vector}
#'   \item{week}{a numeric vector; 0 represents the baseline over the
#'   preceeding 8 weeks}
#'   \item{nSeizures}{a numeric vector; number of seizures}
#'   \item{response}{a factor indicating if the subject had 5 or less seisures.}
#'   }
#' @details
#' Participants in the study were randomized to either progabide or a placebo,
#' as an adjuvant to the standard anti-epileptic chemotherapy. Progabide is an
#' anti-epileptic drug whose primary mechanism of action is to enhance
#' gamma-aminobutyric acid (GABA) content; GABA is the primary inhibitory
#' neurotransmitter in the brain. Prior to receiving treatment, baseline data
#' on the number of epileptic seizures during the preceding 8-week interval were
#' recorded. Counts of epileptic seizures during 2-week intervals before each of
#' four successive post-randomization clinic visits were recorded.
#'
#' @references Thall PF, Vail SC (1990) Some covariance models for longitudinal
#' count data with overdispersion. Biometrics 46:657-671
"epilepsy"



#' Clinical trial comparing two treatments for a respiratory illness
#'
#' The data are from a clinical trial of patients with respiratory illness,
#' where 111 patients from two different clinics were randomized to receive
#' either placebo or an active treatment.
#'
#' @format A data frame with 555 rows and 7 columns:
#' \describe{
#'   \item{id}{a factor with 56 levels}
#'   \item{gender}{a factor with levels F and M}
#'   \item{center}{a factor with levels C1 and C2}
#'   \item{treatment}{a factor with levels placebo and active}
#'   \item{age}{a numeric vector}
#'   \item{visit}{a numeric vector.}
#'   \item{status}{a factor with levels poor and good (response)}
#'   }
#' @details
#' Patients were examined at baseline and at four visits during treatment. At
#' each examination, respiratory status (categorized as 1 = good, 0 = poor) was
#' determined.
#'
#' @references Stokes ME, Davis CS, Koch GG (1995) Categorical Data Analysis
#' using the SAS System. Cary, NC: SAS Institute, Inc.
"respiratory"
