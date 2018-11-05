# this script correct raw expression data by removing principal components

library(argparser)
library(data.table)
library(svd)
library(preprocessCore)
library(limma)

args <- arg_parser('program')

args <- add_argument(args, '-expr',
                     help='expression data file (gene x sample)',
                     default='/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/WholeBlood.txt')
args <- add_argument(args, '-pc',
                     help='number of PCs to remove',
                     default=5)
args <- add_argument(args, '-st',
                     help='column number actual data start at (note: gene names at first col)',
                     default=3)
args <- add_argument(args, '-log',
                     help='log transform',
                     default=TRUE)
args <- add_argument(args, '-quantile',
                     help='quantile transform',
                     default=TRUE)
args <- add_argument(args, '-tpm',
                     help='min given tpm in at least some number of samples',
                     default=1)
args <- add_argument(args, '-sample',
                     help='min tpm in at least given number of samples',
                     default=10)
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/pc_corrected_expr.txt')

argv = parse_args(args)
expr_fn = argv$expr
start_col = argv$st
n_pc = argv$pc
do_log = argv$log
do_quantile = argv$quantile
min_tpm = argv$tpm
min_samples = argv$sample
out_fn = argv$o


### read util
read_df <- function(fn, sep = '\t',  header = T, quote = "", row.names=T, stringsAsFactors = F, check.names = F, lessColInHeader=F, skip = 0){
  if(header==T && lessColInHeader==T){
    header_line = readLines(fn, n = 1)
    headers = strsplit(header_line, split = sep)[[1]]
    skip = skip + 1
  }
  
  data_df = fread(fn, 
                  sep = sep,
                  header = header & !lessColInHeader, 
                  skip = skip,
                  quote = quote,    # not available in old package
                  stringsAsFactors = stringsAsFactors, 
                  check.names = check.names, 
                  data.table = FALSE)
  if (row.names == TRUE){
    rownames(data_df) = data_df[,1]
    data_df = data_df[,-1, drop=F]
  }
  if(header==T && lessColInHeader==T){
    colnames(data_df) = headers
  }
  return(data_df)
}

### read expression
expr_df = read_df(fn = expr_fn, sep = '\t', header = T, quote = "", row.names = T, stringsAsFactors = F, lessColInHeader = F)
start_col_df = NULL
if (start_col>2){
  start_col_df = expr_df[,1:(start_col-2), drop=F]
  expr_df = expr_df[,-(1:(start_col-2))]
}

### filter on TPM
filter_on_tpm_read <- function(expr.df, tpm.df, min.tpm = 0.1, min.samples = 10){
  # expr.df gene x sample dataframe
  
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  tpm.df <- tpm.df[features, samples]
  
  ### filter
  n_samples_w_min_tpm <- apply(tpm.df>min.tpm, 1, sum)
  has.min.samples <- (n_samples_w_min_tpm >= min.samples)
  features.passed <- names(has.min.samples[has.min.samples])
  expr.df <- expr.df[features.passed,]
  
  return(expr.df)
}

expr_df = filter_on_tpm_read(expr.df = expr_df, tpm.df = expr_df, min.tpm = min_tpm, min.samples = min_samples)

### convert to sample x gene matrix
expr_mat = t(expr_df)  # sample x gene

### log
if(do_log == TRUE){
  expr_mat = log2(1+expr_mat)
}

### quantile normalization
qn_expr <- function(expr_df){
  # expr_df: sample x gene dataframe or matrix
  rows = rownames(expr_df)
  cols = colnames(expr_df)
  expr_df = t(normalize.quantiles(t(expr_df)))
  rownames(expr_df) = rows
  colnames(expr_df) = cols
  return(expr_df)
}

if(do_quantile == TRUE){
  expr_mat = qn_expr(expr_mat)
}

### scale
expr_mat = scale(expr_mat)

### pc
expr_svd = propack.svd(expr_mat, neig = n_pc)
pcs = expr_svd$u

### remove pcs
fit = lmFit(t(expr_mat),pcs)
residual_mat = expr_mat - pcs %*% t(fit$coefficients)

### save
residual_df = as.data.frame(t(residual_mat))
if(!is.null(start_col_df))
  residual_df = cbind(start_col_df[rownames(residual_df),,drop=F], residual_df)
write.table(residual_df, file = out_fn, sep = '\t', quote = F, row.names = T, col.names = NA)
