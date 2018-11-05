library(argparser)
library(data.table)
library(igraph)
library(svd)

args <- arg_parser("program");
args <- add_argument(args, '-dir',
                     help='trans-eqtl directory',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/trans_eqtl_cross_chr_mappability_0.8')
args <- add_argument(args, '-tissues',
                     help='tissues -- comma separated',
                     default='Whole_Blood,Muscle_Skeletal,Skin_Sun_Exposed_Lower_leg,Thyroid,Testis')
args <- add_argument(args, '-pfx',
                     help='file name prefix',
                     default='')
args <- add_argument(args, '-sfx',
                     help='file name suffix',
                     default='_cross_chr_map_trans_eqtl_fdr_0.05.txt')
args <- add_argument(args, '-o',
                     help='output file prefix',
                     default="results/combined_trans_eqtl")

argv = parse_args(args)
eqtl_dir = argv$dir
tissues_input = argv$tissues
fn_prefix = argv$pfx
fn_suffix = argv$sfx
output_prefix = argv$o


### directory must exist
stopifnot(dir.exists(eqtl_dir))

### parse tissues
parts = strsplit(tissues_input, split = ',')[[1]]
tissues = parts[sapply(parts, nchar)>0]

### read and combine files
eqtl_dfs = lapply(tissues, function(tissue){
  fn = sprintf("%s/%s%s%s", eqtl_dir, fn_prefix, tissue, fn_suffix)
  stopifnot(file.exists(fn))
  eqtl_df = read.table(fn, header = T, sep = '\t', quote = "", comment.char = "", stringsAsFactors = F, colClasses = 'character')
  eqtl_df$tissue = tissue
  return(eqtl_df)
})
combined_eqtl_df = do.call(rbind, eqtl_dfs)
uniq_combined_eqtl_df = aggregate( . ~ snps + gene, data = combined_eqtl_df,  FUN=function(vals) paste0(vals, collapse = ','))

### save
combined_fn = sprintf("%s.all.txt", output_prefix)
uniqu_combined_fn = sprintf("%s.all.unique.txt", output_prefix)
write.table(combined_eqtl_df, file = combined_fn, sep = '\t', row.names = F, col.names = T, quote = F)
write.table(uniq_combined_eqtl_df, file = uniqu_combined_fn, sep = '\t', row.names = F, col.names = T, quote = F)
