test_that("Table 1 of Touloumis (2026a) is reproduced (shoulder-tip pain)", {
  fits <- lapply(
    c(gee = "gee", rbr = "brgee-robust", rbc = "bcgee-robust"),
    function(method) {
      geewa_binary(
        pain ~ suction + age + female + last_day,
        id = id,
        repeated = time,
        data = test_data$cholecystectomy,
        link = "logit",
        orstr = "unstructured",
        method = method
      )
    }
  )
  expect_published(min(fits$gee$alpha), 3.99, 2)
  expect_published(max(fits$gee$alpha), 135.00, 2)
  expect_published_table(
    fits$gee,
    c(-2.1898, 1.8602, 0.0317, 0.3478, 0.9513),
    c(1.0698, 0.5794, 0.0165, 0.5953, 0.3894),
    c(0.0407, 0.0013, 0.0548, 0.5591, 0.0146)
  )
  expect_published_table(
    fits$rbr,
    c(-2.0360, 1.7602, 0.0295, 0.3065, 0.9029),
    c(1.0233, 0.5566, 0.0159, 0.5773, 0.3770),
    c(0.0466, 0.0016, 0.0626, 0.5955, 0.0166)
  )
  expect_published_table(
    fits$rbc,
    c(-2.0118, 1.7583, 0.0291, 0.2925, 0.9079),
    c(1.0188, 0.5549, 0.0158, 0.5760, 0.3766),
    c(0.0483, 0.0015, 0.0654, 0.6116, 0.0159)
  )
})

test_that("Table 1 of Touloumis (2026b) is reproduced (respiratory trial)", {
  fits <- lapply(
    c(
      independence = "independence",
      exchangeable = "exchangeable",
      unstructured = "unstructured"
    ),
    function(orstr) {
      geewa_binary(
        status ~ active + male + visit + age + baseline + active:male +
          visit:age,
        id = id,
        repeated = visit,
        data = test_data$respiratory2,
        link = "probit",
        orstr = orstr,
        method = "pgee-jeffreys"
      )
    }
  )
  for (fit in fits) {
    expect_lt(summary(fit)$coefficients["active", "Pr(>|z|)"], 1e-4)
  }
  expect_published(min(fits$unstructured$alpha), 3.56, 2)
  expect_published(max(fits$unstructured$alpha), 19.04, 2)
  expect_published_table(
    fits$independence,
    c(-1.1848, 2.3727, 0.2927, 0.3542, 0.0203, 0.6650, -1.7077, -0.0106),
    c(0.7797, 0.5196, 0.4345, 0.2204, 0.0162, 0.3175, 0.6293, 0.0053),
    c(0.1286, NA, 0.5005, 0.1080, 0.2103, 0.0362, 0.0067, 0.0449)
  )
  expect_published_table(
    fits$exchangeable,
    c(-1.3384, 2.3341, 0.2795, 0.3994, 0.0250, 0.6604, -1.6485, -0.0119),
    c(0.7559, 0.5096, 0.4458, 0.2044, 0.0155, 0.3194, 0.6228, 0.0049),
    c(0.0766, NA, 0.5307, 0.0507, 0.1082, 0.0387, 0.0081, 0.0145)
  )
  expect_published_table(
    fits$unstructured,
    c(-1.2292, 2.3229, 0.3130, 0.3509, 0.0195, 0.6976, -1.6341, -0.0100),
    c(0.7164, 0.5129, 0.4383, 0.1957, 0.0146, 0.3136, 0.6221, 0.0046),
    c(0.0862, NA, 0.4751, 0.0730, 0.1803, 0.0261, 0.0086, 0.0295)
  )
})
