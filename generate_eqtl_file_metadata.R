library(argparser)
library(data.table)
library(parallel)

args <- arg_parser('program')
args <- add_argument(args, '-pfx', 
                     help='file prefix for eqtl results file',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls/Whole_Blood/meqtl_Whole_Blood')
args <- add_argument(args, '-core',
                     help='number of cores',
                     default=1)
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/metadata.txt')

argv <- parse_args(args)
file_prefix = argv$pfx
n_cores = argv$core
out_fn <- argv$o


### get file names
file_dir = dirname(file_prefix)
pattern = paste0(basename(file_prefix), '.*.RData')
files = list.files(path = file_dir, pattern = pattern, full.names = T)

### for each file, read and find chr, left snp pos, right snp pos
meta_data_list = mclapply(files, mc.cores = n_cores, FUN = function(fn){
  print(fn)
  load(fn)
  snps = levels(meqtl_df$snps)
  snps_pos = data.frame(t(sapply(snps, function(s) strsplit(s, '_')[[1]][1:2])), stringsAsFactors = F)
  stopifnot(length(unique(snps_pos[,1]))==1)
  snps_pos[,2] = as.integer(snps_pos[,2])
  metadata = c(file=basename(fn), chr=as.character(snps_pos[1,1]), start_pos=min(snps_pos[,2]), end_pos=max(snps_pos[,2]))
  rm(meqtl_df)
  return(metadata)
})

meta_data_df = do.call(rbind,  meta_data_list)

### write output
write.table(meta_data_df, file = out_fn, sep = '\t', quote = F, row.names = F, col.names = T)
