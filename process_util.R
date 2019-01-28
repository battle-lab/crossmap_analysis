library('edgeR')

to_inv_normal <- function(expr_mat){
  # expr_mat : gene x sample matrix
  transformed_data <- apply(expr_mat, 1, function(x){
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  })
  transformed_data = t(transformed_data)
  return(transformed_data)
}

tmm_normalize <- function(expr_mat){
  # expr_mat: gene x sample matrix
  tmm_factors = calcNormFactors(expr_mat, method = "TMM")
  expr_mat = t(t(expr_mat)/tmm_factors)
  return(expr_mat)
}
