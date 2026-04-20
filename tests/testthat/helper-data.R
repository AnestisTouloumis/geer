local_data <- function() {
  data("epilepsy", package = "geer", envir = environment())
  data("cerebrovascular", package = "geer", envir = environment())
  data("respiratory", package = "geer", envir = environment())
  data("leprosy", package = "geer", envir = environment())
  data("rinse", package = "geer", envir = environment())
  list(
    epilepsy = epilepsy,
    cerebrovascular = cerebrovascular,
    respiratory = respiratory,
    respiratory2 = respiratory[respiratory$center == "C2", , drop = FALSE],
    leprosy = leprosy,
    rinse = rinse
  )
}
test_data <- local_data()
