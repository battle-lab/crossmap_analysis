library(argparser)
library(data.table)
library(parallel)

args <- arg_parser('program')
args <- add_argument(args, "-files", 
                     help="shell script file patterns for matrix-eQTL results", 
                     default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls/Whole_Blood/meqtl_Whole_Blood*.RData")
args <- add_argument(args, '-annot',
                     help='gene annotation file',
                     default='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, '-d',
                     help='distance between gene and snp',
                     default=1e6)
args <- add_argument(args, '-core',
                     help='number of cores.',
                     default=1)
args <- add_argument(args, '-o',
                     help='output file.',
                     default='results/best_cis_snp_per_gene.txt')


argv <- parse_args(args)
meqtl_files_pattern = argv$files
gene_annot_fn = argv$annot
d_gene_snp = argv$d
n_cores = argv$core
out_fn <- argv$o

meqtl_files = system(sprintf('ls %s', meqtl_files_pattern), intern = T)
stopifnot(is.null(attr(meqtl_files, 'status')))
stopifnot(length(meqtl_files)>0)


### read gene annotation file and compute tss
gene_annot_df = read.table(gene_annot_fn, header = T, sep = '\t', row.names = 1, quote = "", comment.char = "", stringsAsFactors = F)
dim(gene_annot_df)
head(gene_annot_df)

tss =  as.integer(apply(gene_annot_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))
gene_annot_df$tss = as.integer(tss)

### function to load matrix-eqtl results
get_meqtl_results <- function(fn){
  # load and insert the file
  stopifnot(file.exists(fn))
  load(fn)
  ## add snp chr and snp pos
  meqtl_snps_levels = levels(meqtl_df$snps)
  meqtl_snps_levels_chr = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) gsub(pattern = 'chr', replacement = '', x = parts[1], ignore.case = T) )
  meqtl_df$snp_chr = meqtl_snps_levels_chr[as.integer(meqtl_df$snps)]
  
  meqtl_snps_levels_pos = as.integer(sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) parts[2]))
  meqtl_df$snp_pos = meqtl_snps_levels_pos[as.integer(meqtl_df$snps)]
  
  ## add gene chr and tss
  meqtl_gene_levels = levels(meqtl_df$gene)
  meqtl_gene_levels_chr = sapply(gene_annot_df[meqtl_gene_levels,'chr'], function(x) gsub(pattern = 'chr', replacement = '', x = x, ignore.case = T))
  meqtl_df$gene_chr = meqtl_gene_levels_chr[as.integer(meqtl_df$gene)]
  
  meqtl_gene_levels_tss = gene_annot_df[meqtl_gene_levels,'tss']
  meqtl_df$gene_tss = meqtl_gene_levels_tss[as.integer(meqtl_df$gene)]
  
  return(meqtl_df)
}

best_snp_per_gene_in_dataframe <- function(res_meqtl_df){
  candidates = res_meqtl_df[res_meqtl_df$snp_chr == res_meqtl_df$gene_chr, ]
  candidates$d = abs(candidates$snp_pos - candidates$gene_tss)
  candidates = candidates[candidates$d <= d_gene_snp, ]
  candidates = candidates[with(candidates, order(gene, pvalue, d)), ]
  best_snp_per_file = candidates[!duplicated(candidates$gene),]
  rm(candidates)
  return(best_snp_per_file)
}

best_snp_per_gene_per_file <- mclapply(X = meqtl_files, mc.cores = n_cores, FUN = function(fn){
  print(fn)
  res_meqtl_df = get_meqtl_results(fn)
  best = best_snp_per_gene_in_dataframe(res_meqtl_df)
  rm(res_meqtl_df)
  gc(reset = T)
  return(best)
})

combined_best_snps_per_file_df <- do.call(rbind, best_snp_per_gene_per_file)
combined_best_snp_per_gene <- best_snp_per_gene_in_dataframe(combined_best_snps_per_file_df)

write.table(combined_best_snp_per_gene, file = out_fn, sep = '\t', quote = F, row.names = F, col.names = T)
