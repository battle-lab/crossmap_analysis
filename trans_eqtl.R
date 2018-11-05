library(argparser)
library(MatrixEQTL)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program to run matrix-eqtl");
args <- add_argument(args, "-expr", help="expression file", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Muscle_Skeletal.v7.normalized_expression.txt")
args <- add_argument(args, "-snp", help="snp file", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process/chr1_maf_0.05_repeat_masked.012.txt")
args <- add_argument(args, "-cov", help="expression file", default="/work-zfs/abattle4/lab_data/GTEx_v7_cisEQTL/GTEx_Analysis_v7_eQTL_covariates/Muscle_Skeletal.v7.covariates.txt")
args <- add_argument(args, "-max_snp", help="max snp per matrix-eqtl call", default=1000)
args <- add_argument(args, "-o", help="out prefix (RData)", default="results/eqtl")

argv = parse_args(args)
expr_fn = argv$expr
snp_fn = argv$snp
cov_fn = argv$cov
max_snp_per_meqtl = argv$max_snp
out_prefix = argv$o

### read data
expr_mat = as.matrix(read_df(expr_fn))

snp_mat = as.matrix(read_df(snp_fn))
snp_mat[snp_mat<0] = NA

cov_df = as.data.frame(t(read_df(cov_fn)))
n_uniq = sapply(cov_df, function(x)length(unique(x)))
for(covariate in names(n_uniq[n_uniq<=5])){
  cov_df[,covariate] = as.factor(cov_df[,covariate])
}
#remove constant covariates
if(any(n_uniq<=1)){ 
  const_covariates = names(n_uniq[n_uniq<=1])
  cov_df = cov_df[,-which(colnames(cov_df) %in% const_covariates), drop=F]
}

common_samples = intersect(colnames(expr_mat), colnames(snp_mat))
common_samples = intersect(common_samples, rownames(cov_df))

expr_mat = expr_mat[,common_samples]
snp_mat = snp_mat[,common_samples]
cov_mat = t(model.matrix( ~ ., data=cov_df[common_samples,]))
cov_mat = cov_mat[-which(rownames(cov_mat)=='(Intercept)'),]

### run matrix-eQTL : part by part
meqtl_expr_slice = SlicedData$new(expr_mat)
meqtl_cov_slice = SlicedData$new(cov_mat)

snp_splits = sort(unique(c(seq(max_snp_per_meqtl, nrow(snp_mat), max_snp_per_meqtl), nrow(snp_mat))))
snp_splits = snp_splits[snp_splits<=nrow(snp_mat)]

for(ss in 1:length(snp_splits)){
  print(sprintf("calling matrix-eqtl: part %d of %d", ss, length(snp_splits) ))
  start_idx = ifelse(ss==1, 1, snp_splits[ss-1]+1)
  end_idx = snp_splits[ss]
  meqtl_snp_slice = SlicedData$new(snp_mat[start_idx:end_idx,,drop=F])
  
  me = Matrix_eQTL_engine(
    snps = meqtl_snp_slice,
    gene = meqtl_expr_slice,
    cvrt = meqtl_cov_slice,
    output_file_name = NULL,
    pvOutputThreshold = 1,
    useModel = modelLINEAR, 
    verbose = FALSE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  meqtl_df = me$all$eqtls
  eqtl_fn = sprintf("%s_part_%d.RData", out_prefix, ss)
  save(meqtl_df, file = eqtl_fn)
  
  rm(meqtl_snp_slice, me, meqtl_df)
  gc(reset = T)
  
}

