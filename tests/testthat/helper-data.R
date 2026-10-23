load_geer_data <- function(name) {
  env <- new.env(parent = emptyenv())
  data(list = name, package = "geer", envir = env)
  get(name, envir = env, inherits = FALSE)
}

build_agresti_data <- function() {
  counts <- rbind(
    c(16, 13, 9, 3, 14, 4, 15, 6),
    c(31, 0, 6, 0, 22, 2, 9, 0),
    c(2, 2, 8, 9, 9, 15, 27, 28),
    c(7, 2, 5, 2, 31, 5, 32, 6)
  )
  patterns <- c("NNN", "NNA", "NAN", "NAA", "ANN", "ANA", "AAN", "AAA")
  groups <- data.frame(diagnosis = c(0, 0, 1, 1), treatment = c(0, 1, 0, 1))
  subjects <- do.call(rbind, lapply(seq_len(nrow(groups)), function(g) {
    data.frame(
      diagnosis = groups$diagnosis[g],
      treatment = groups$treatment[g],
      pattern = rep(patterns, counts[g, ])
    )
  }))
  data.frame(
    id = rep(seq_len(nrow(subjects)), each = 3),
    time = rep(0:2, nrow(subjects)),
    diagnosis = rep(subjects$diagnosis, each = 3),
    treatment = rep(subjects$treatment, each = 3),
    abnormal = as.numeric(unlist(strsplit(subjects$pattern, "")) == "A")
  )
}

build_cholecystectomy_data <- function() {
  data <- load_geer_data("cholecystectomy")
  data$suction <- as.numeric(data$treatment == "active")
  data$female <- as.numeric(data$gender == "female")
  data$last_day <- as.numeric(data$time >= 5)
  data
}

build_respiratory2_data <- function(respiratory) {
  data <- respiratory[respiratory$center == "C2", , drop = FALSE]
  data$active <- as.numeric(data$treatment == "active")
  data$male <- as.numeric(data$gender == "male")
  data
}

test_data <- local({
  respiratory <- load_geer_data("respiratory")
  list(
    agresti = build_agresti_data(),
    cerebrovascular = load_geer_data("cerebrovascular"),
    cholecystectomy = build_cholecystectomy_data(),
    epilepsy = load_geer_data("epilepsy"),
    leprosy = load_geer_data("leprosy"),
    respiratory = respiratory,
    respiratory2 = build_respiratory2_data(respiratory),
    rinse = load_geer_data("rinse")
  )
})
