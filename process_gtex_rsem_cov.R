suppressMessages(library(peer))
suppressMessages(library('argparser'))
suppressMessages(source('io_util.R'))

args <- arg_parser('program')
args <- add_argument(args, "-expr", 
                     help="expression file (gene x sample)", 
                     default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/processed_expression/Whole_Blood.v7.normalized_expression.txt")
args <- add_argument(args, "-prev", 
                     help="previous covariate file (cov x sample)", 
                     default="/work-zfs/abattle4/lab_data/GTEx_v7_cisEQTL/GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt")
args <- add_argument(args, '-o',
                     help='output covariate file.',
                     default='results/cov.txt')


argv <- parse_args(args)
expr_fn = argv$expr
prev_fn = argv$prev
out_fn <- argv$o

prev_cov_variables = c('C1', 'C2', 'C3', 'sex', 'platform')

### read data
expr_df = read_df(expr_fn)
prev_df = read_df(prev_fn)

expr_mat = t(expr_df)
prev_cov_mat = t(prev_df[intersect(rownames(prev_df), prev_cov_variables), colnames(expr_df)])

### infer peer factors
n_factor = nrow(prev_df) - ncol(prev_cov_mat)
model = PEER()
PEER_setPhenoMean(model,expr_mat)
PEER_setCovariates(model, prev_cov_mat)
PEER_setNk(model,n_factor)
PEER_setNmax_iterations(model, 1000)
PEER_update(model)
factors = PEER_getX(model)

### organize factor matrix
rownames(factors) = rownames(expr_mat)
colnames(factors) = c(colnames(prev_cov_mat), paste0('InferredCov', 1:n_factor))
factors = t(factors) # cov x sample

### save
write_df(factors, file = out_fn)

