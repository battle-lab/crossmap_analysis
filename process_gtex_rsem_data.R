suppressMessages(source('io_util.R'))
suppressMessages(source('process_util.R'))
suppressMessages(library('argparser'))


args <- arg_parser('program')
args <- add_argument(args, "-tpm", 
                     help="tpm file (gene x sample)", 
                     default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/raw/Whole_Blood.tpm.txt")
args <- add_argument(args, "-prev", 
                     help="previous processed file (gene x sample)", 
                     default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt")
args <- add_argument(args, '-o',
                     help='output normalized file.',
                     default='results/normalized.txt')


argv <- parse_args(args)
tpm_fn = argv$tpm
prev_fn = argv$prev
out_fn <- argv$o

### read data
tpm_df = read_df(tpm_fn)
prev_df = read_df(prev_fn)

### take samples and genes from prev data
indiv_ids = sapply(colnames(tpm_df), function(s) paste0(strsplit(s,split='-')[[1]][1:2], collapse = '-') )
colnames(tpm_df) = indiv_ids
# new data must have all rows and columns of prev data
stopifnot(length(setdiff(rownames(prev_df), rownames(tpm_df))) == 0)
stopifnot(length(setdiff(colnames(prev_df), colnames(tpm_df))) == 0)
tpm_df = tpm_df[rownames(prev_df),colnames(prev_df)]

### normalize using TMM
tpm_df = tmm_normalize(tpm_df)

### inverse normal
tpm_df =  to_inv_normal(tpm_df)

### save
write_df(tpm_df, file = out_fn)

