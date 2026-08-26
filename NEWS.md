# geer 0.1.1

## New features

* Added optional integration with the `marginaleffects` package. `geer` model
  objects can now be used with `marginaleffects` workflows for predictions,
  comparisons, slopes, and their averaged counterparts. The integration uses
  delayed S3 registration, so `marginaleffects` remains an optional dependency,
  and includes model-data recovery through `insight`.

* Expanded `geecriteria()` with additional criteria for working association
  structure and model selection: `QICHH`, `QICC`, `EQIC`, `GHYC`, `PAC`,
  `AGPC`, and `SGPC`. These complement the existing `QIC`, `CIC`, `RJC`,
  `QICu`, `GESSC`, and `GPC` criteria.

* Added `runs_test()` for the Wald-Wolfowitz nonparametric runs test of GEE
  residual signs described by Chang (2000) and Hardin and Hilbe (2013). The
  test supports Pearson, deviance, and working residuals and can assess the
  natural cluster/repeated order, fitted-value order, or covariate-based
  orderings.

* Changed the default covariance estimator in `geecriteria()` to `"robust"` so
  the classical forms of covariance-based GEE criteria are returned by default.
  Other covariance estimators remain available through `cov_type`.
