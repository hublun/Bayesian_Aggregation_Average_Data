`.sourceCpp_1_DLLInfo` <- dyn.load('/home/dulunche/Bayesian_Aggregation_Average_Data/rcpp-cache/sourceCpp-x86_64-pc-linux-gnu-1.0.12/sourcecpp_2217dd27699/sourceCpp_2.so')

valogit <- Rcpp:::sourceCppFunction(function(va, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_valogit')
inv_valogit <- Rcpp:::sourceCppFunction(function(tva, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_inv_valogit')
evaluate_model <- Rcpp:::sourceCppFunction(function(x, DRUG, Lalpha_0, Lalpha_s, lkappa, lEmax, delta, eta_Lalpha_0, eta_lkappa, sigma_Lalpha_0, sigma_lkappa, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_evaluate_model')
pretty_print <- Rcpp:::sourceCppFunction(function(x, pstream__ = 0L) {}, TRUE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_pretty_print')
colMeans <- Rcpp:::sourceCppFunction(function(y, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_colMeans')
robust_cov <- Rcpp:::sourceCppFunction(function(y, ym, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_robust_cov')
mvn_approx_lpdf <- Rcpp:::sourceCppFunction(function(y_prime_bar, x_prime, J_prime, DRUG_prime, delta, Lalpha_0, Lalpha_s, lkappa, lEmax, sigma_Lalpha_0, sigma_lkappa, sigma_y, xi, pstream__ = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mvn_approx_lpdf')

rm(`.sourceCpp_1_DLLInfo`)
