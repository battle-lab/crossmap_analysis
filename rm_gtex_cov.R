# this script correct raw expression data by removing principal components

library(argparser)
library(limma)

args <- arg_parser('program')

args <- add_argument(args, '-expr',
                     help='expression data file (gene x sample)',
                     default='/work-zfs/abattle4/lab_data/GTEx_v7_cisEQTL/GTEx_Analysis_v7_eQTL_expression_matrices/Whole_Blood.v7.normalized_expression.bed.gz')
args <- add_argument(args, '-cov',
                     help='covariate daa file (cov x sample)',
                     default='/work-zfs/abattle4/lab_data/GTEx_v7_cisEQTL/GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt')
args <- add_argument(args, '-id',
                     help='column number for gene id in expression file',
                     default=4)
args <- add_argument(args, '-st',
                     help='column number actual data start at in expression file',
                     default=5)
args <- add_argument(args, '-o',
                     help='corrected output file (gene x sample)',
                     default='results/corrected.txt')

argv = parse_args(args)
expr_fn = argv$expr
cov_fn = argv$cov
id_col = argv$id
start_col = argv$st
out_fn = argv$o

stopifnot(start_col>=2)
stopifnot(id_col>=1)
stopifnot(id_col<start_col)

### read data
expr_df = read.table(file = expr_fn, sep = '\t', header = T, comment.char = "", quote = "", stringsAsFactors = F, check.names = F)
rownames(expr_df) = expr_df[,id_col]
expr_val_df = expr_df[,start_col:ncol(expr_df)]

start_col_df = expr_df[,1:(start_col-1), drop=F]


cov_df = read.table(file = cov_fn, sep = '\t', header = T, comment.char = "", quote = "", stringsAsFactors = F, row.names = 1, check.names = F)
stopifnot(length(intersect(colnames(expr_val_df), colnames(cov_df))) == ncol(expr_val_df))
cov_df = cov_df[,colnames(expr_val_df)]

cov_df = as.data.frame(t(cov_df))            # sample x cov
n_uniq = sapply(cov_df, function(x)length(unique(x)))
for(covariate in names(n_uniq[n_uniq<=5])){
  cov_df[,covariate] = as.factor(cov_df[,covariate])
}
#remove constant covariates
if(any(n_uniq<=1)){ 
  const_covariates = names(n_uniq[n_uniq<=1])
  cov_df = cov_df[,-which(colnames(cov_df) %in% const_covariates), drop=F]
}


### remove covariates
expr_mat = as.matrix(expr_val_df)  # gene x sample
design_mat = model.matrix( ~ ., data=cov_df)
fit = lmFit(expr_mat, design_mat)
residual_mat = expr_mat - fit$coefficients %*% t(design_mat)

### save
residual_df = as.data.frame(residual_mat)
residual_df = cbind(start_col_df[rownames(residual_df),,drop=F], residual_df)
write.table(residual_df, file = out_fn, sep = '\t', quote = F, row.names = F, col.names = T)
