load_geer_data <- function(name) {
  env <- new.env(parent = emptyenv())
  data(list = name, package = "geer", envir = env)
  get(name, envir = env, inherits = FALSE)
}

test_data <- local({
  respiratory <- load_geer_data("respiratory")
  list(
    cerebrovascular = load_geer_data("cerebrovascular"),
    epilepsy = load_geer_data("epilepsy"),
    leprosy = load_geer_data("leprosy"),
    respiratory = respiratory,
    respiratory2 = respiratory[respiratory$center == "C2", , drop = FALSE],
    rinse = load_geer_data("rinse")
  )
})
