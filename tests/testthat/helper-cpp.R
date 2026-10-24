get_cached_cpp_utils <- local({
  functions <- new.env(parent = baseenv())
  function() {
    if (!exists("compute_cpp_utils", envir = functions, inherits = FALSE)) {
      source <- testthat::test_path("..", "..", "src", "utils.cpp")
      if (!file.exists(source)) {
        testthat::skip("the C++ utility source is unavailable")
      }
      source <- normalizePath(source, winslash = "/")
      Rcpp::cppFunction(
        depends = "RcppArmadillo",
        plugins = "cpp17",
        includes = paste0('#include "', source, '"'),
        code = paste(
          "arma::mat compute_cpp_utils(arma::mat x, const arma::vec& b,",
          "                            const int route) {",
          "  if (route == 2) {",
          '    symmetrize_if_close(x, "compute_cpp_utils");',
          "    return x;",
          "  }",
          '  const CholOrLuFactor factor(x, "compute_cpp_utils", route == 0);',
          "  return factor.solve(b);",
          "}",
          sep = "\n"
        ),
        env = functions
      )
    }
    get("compute_cpp_utils", envir = functions, inherits = FALSE)
  }
})
