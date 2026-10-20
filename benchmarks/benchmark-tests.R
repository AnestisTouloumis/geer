## geer
devtools::document()
devtools::test()

## gee -- all pass
testthat::test_file("benchmarks/benchmark-gee.R")

## geeM -- all pass
testthat::test_file("benchmarks/benchmark-geeM.R")

## BCgee -- 1 skip -- offset
testthat::test_file("benchmarks/benchmark-BCgee.R")

## geessbin -- 9 skip -- centered Morel et al.
testthat::test_file("benchmarks/benchmark-geessbin.R")

## geepack -- 13 skip -- not comparable (8) & QICC (4) & error code (1)
testthat::test_file("benchmarks/benchmark-geepack.R")

## glmtoolbox -- 24 skip -- balanced designs (18) & not comparable (6)
testthat::test_file("benchmarks/benchmark-glmtoolbox.R")
