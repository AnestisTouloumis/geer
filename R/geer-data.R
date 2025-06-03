#' Cerebrovascular Deficiency Crossover Trial
#'
#' The dataset consists of safety data from a crossover trial on cerebrovascular deficiency.
#'
#' @format A data frame with 134 rows and 4 columns:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{period}{Period identifier.}
#'   \item{ecg}{Response indicating whether an electrocardiogram was
#'   abnormal or normal.}
#'   \item{treatment}{Treatment group variable with two levels: placebo and active.}
#'   }
#' @details
#' The response variable is not a trial endpoint but rather a potential side
#' effect, indicating whether an electrocardiogram (ECG) was abnormal or normal.
#' This two-period crossover trial compares the effects of active drug to placebo.
#' The 67 patients were randomly allocated to the two treatment sequences, with
#' 34 patients receiving placebo followed by active treatment, and 33 patients
#' receiving active treatment followed by placebo.
#'
#' @references Jones, B. and Kenward, M.G. (1989) \emph{Design and Analysis of Cross-over
#' Trials.} London: Chapman and Hall/CRC Press.
#' data("cerebrovascular")
#' str("cerebrovascular")
"cerebrovascular"


#' Epilepsy Clinical Trial
#'
#' The data are from a placebo-controlled clinical trial on epilepsy.
#'
#' @format A data frame with 295 rows and 6 columns:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{week}{Week identifier: 0 represents the baseline over the
#'   preceeding 8 weeks}
#'   \item{y}{Response coded as (1) if the subject had 5 or less seizures and as (0) otherwise.}
#'   \item{treatment}{Treatment group variable with two levels: placebo and progabide.}
#'   \item{age}{Age recorded at baseline.}
#'   \item{nSeizures}{The number of seizures.}
#'   }
#' @details
#' A total of 57 patients with partial seizures were enrolled in a randomized
#' clinical trial. Patients were allocated to either the anti-epileptic drug
#' progabide or to placebo. Prior to receiving treatment,
#' the number of epileptic seizures during the preceding 8-week interval were
#' recorded. After receiving treatment, the number of epileptic seizures were
#' recorded at each of four 2-week intervals clinic visits.
#'
#' @references Thall, P.F. and Vail, S.C. (1990) Some covariance models for longitudinal
#' count data with overdispersion. \emph{Biometrics}, \strong{46}, 657-671.
#' @examples
#' data("epilepsy")
#' str("epilepsy")
"epilepsy"



#' Respiratory Clinical Trial
#'
#' The data are from a clinical trial of patients with respiratory illness.
#'
#' @format A data frame with 444 rows and 7 columns:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{y}{Respiratory status variable with two levels: good and poor.}
#'   \item{baseline}{Respiratory status variable at baseline with two levels:
#'    good and poor.}
#'   \item{treatment}{Treatment group variable with two levels: placebo and active.}
#'   \item{visit}{Visit identifier.}
#'   \item{gender}{The gender of the patient with two levels: male and female.}
#'   \item{age}{The age in years recorded at baseline.}
#'   \item{center}{Center identifier.}
#'   }
#' @details
#' A total of 111 patients from two different centers were enrolled in a randomized
#' clinical trial. Prior to receiving treatment, patients were examined at baseline.
#' After receiving treatment (placebo or active), patients were examined at four
#' visits during. The respiratory status was determined at the baseline and at
#' each of the four follow-up visits.
#'
#' @references Stokes, M.E., Davis, C.S. and Koch, G.G. (1995) \emph{Categorical Data Analysis
#' using the SAS System.} Cary, NC: SAS Institute, Inc.
#' @examples
#' data("respiratory")
#' str("respiratory")
"respiratory"



#' Obesity Study
#'
#' The data are from the Missouri Adolescent Female Twin Study.
#'
#' @format A data frame with 195 rows and 7 columns:
#' \describe{
#'   \item{fid}{Family identifier.}
#'   \item{sid}{Individual identifier.}
#'   \item{obesity}{Obesity status with three levels: Lean, Obese and
#'   Overweight.}
#'   \item{age}{Age in years.}
#'   \item{zygosity}{Zygocity identifier with two levels: monozygotic (MZ) or
#'   dizygotic (DZ).}
#'   \item{ancestry}{Ancestry identifier with two levels: European (EA) or
#'   African (AA).}
#'   \item{bacteroides}{the abundance of bacteroides.}
#'   }
#' @details
#' The data were collected from 54 families with adult female twin pairs born
#' in Missouri. The obesity category for each individual was recorded as Obese
#' (BMI > 29), lean (BMI < 25) or overweight (25 < BMI < 30). One of the aims of
#'  the study was to understand the relationship between obesity and gut microbiome.
#'  Individual RNA samples were obtained at the baseline and at most at one
#'  follow up time and the abundance of several microbiomes was measured
#'  using RNA sequencing. This dataset is a subset of the original dataset that
#'  includes only the abundance of Bacteroides. In addition, the age, the
#'  zygocity and the ancestry of each twin were recorded.
#'
#' @references Turnbaugh et al. (2009) A core gut microbiome in obese and lean twins. \emph{Nature} , \strong{457}, 480–484.
#' @examples
#' data("obesity")
#' str("obesity")
"obesity"


#' Leprosy Clinical Trial
#'
#' The data are from a clinical trial of patients with leprosy.
#'
#' @format A data frame with 80 rows and 4 columns:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{period}{Period identifier.}
#'   \item{time}{Indicator of the post-treatment measurement.}
#'   \item{bacilli}{The number of leprosy bacilli at six sites of the body.}
#'   \item{treatment}{Treatment group variable with three levels: A, B and C.}
#'   }
#' @details
#' A total of 30 patients were enrolled in a randomized clinical trial in
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
#' @examples
#' data("leprosy")
#' str("leprosy")
"leprosy"
