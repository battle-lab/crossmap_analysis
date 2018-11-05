### this script collects eqtl status in different gtex tissues

library(argparser)
source('io_util.R')
source('arg_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="file with all eqtls with crossmapping info", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt")
args <- add_argument(args, "-tis", help="tissues separated by comma", default="Whole_Blood,Muscle_Skeletal,Thyroid,Skin_Sun_Exposed_Lower_leg,Testis")
args <- add_argument(args, "-meqtl", help="matrix eqtl results file format, replace #TISSUE#, #CHR#, #PART#", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls/#TISSUE#/#FILE#")
args <- add_argument(args, "-meqtlmeta", help="matrix eqtl metadata file format, replace #TISSUE#", default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls_metadata/meqtl_meta_#TISSUE#.txt')
args <- add_argument(args, "-o", help="output file", default="results/combined_eqtl_stats.txt")

argv = parse_args(args)
eqtl_fn = argv$eqtl
tissue_input = argv$tis
matrix_eqtl_results_file_format = argv$meqtl
matrix_eqtl_metadata_file_format = argv$meqtlmeta
out_fn = argv$o

tissues = parse_delimitted_param(tissue_input, delim = ',', rm.empty = T)

### read inputs
eqtl_df = read_df(eqtl_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")
  
tissue_stats_per_tissue = lapply(tissues, function(tissue){
  print(sprintf('gathering stats from %s ...', tissue))
  ### read meqtl meta data
  meqtl_meta_fn = gsub(pattern = '#TISSUE#', replacement = tissue, x = matrix_eqtl_metadata_file_format)
  meqtl_meta_df = read_df(meqtl_meta_fn, row.names = F, header = T)
  meqtl_meta_df$chr = paste0('chr', meqtl_meta_df$chr)
  
  ### get matrix-eqtl results file name for each eqtl
  get_meqtl_results_fn <- function(snp_chr, snp_pos, meqtl_meta_df){
    return(meqtl_meta_df[meqtl_meta_df$chr == snp_chr & meqtl_meta_df$start_pos <= snp_pos & meqtl_meta_df$end_pos >= snp_pos, 'file'])
  }
  eqtl_df_tissue = eqtl_df
  eqtl_df_tissue$meqtl_fn = mapply(function(snp_chr, snp_pos){get_meqtl_results_fn(snp_chr, snp_pos, meqtl_meta_df)}, eqtl_df$snps_chr, eqtl_df$snps_pos, SIMPLIFY = T)
  
  ### read matrix-eqtl results
  matrix_eqtl_results_file_format_tissue =gsub(pattern = '#TISSUE#', replacement = tissue, x = matrix_eqtl_results_file_format)
  tissue_stats_per_file = lapply(sort(unique(eqtl_df_tissue$meqtl_fn)), function(fn){
    print(sprintf('gathering stats from %s: %s', tissue, fn))
    meqtl_fn = gsub(pattern = '#FILE#', replacement = fn, x = matrix_eqtl_results_file_format_tissue)
    stopifnot(file.exists(meqtl_fn))
    load(meqtl_fn)
    #meqtl_df$snps = as.character(meqtl_df$snps)
    #meqtl_df$gene = as.character(meqtl_df$gene)
    x1 = eqtl_df_tissue[eqtl_df_tissue$meqtl_fn == fn, c('snps', 'gene')]
    x2 = meqtl_df[meqtl_df$gene %in% unique(x1$gene), ]
    x2 = x2[x2$snps %in% unique(x1$snps), ]
    return(merge(x1, x2, by=c('snps', 'gene'), all.x=T, all.y=F))
  })
  tissue_stats_df = do.call(rbind,tissue_stats_per_file)
  return(tissue_stats_df)
})
names(tissue_stats_per_tissue) = tissues

### combine stats per tissue
combined_tissue_stats = tissue_stats_per_tissue[[1]]
for(idx in 1:length(tissue_stats_per_tissue)){
  combined_tissue_stats = merge(combined_tissue_stats, tissue_stats_per_tissue[[idx]], by=c('snps','gene'), suffixes=c('',paste0('_', tissues[idx]) ))
}
combined_tissue_stats = combined_tissue_stats[,-which(!colnames(tissue_stats_per_tissue[[1]]) %in% c('snps','gene'))]

### merge with given annotated eqtls
eqtl_w_tissue_stats = merge(eqtl_df, combined_tissue_stats, by = c('snps','gene'), all.x=T)

### save
write_df(eqtl_w_tissue_stats, file = out_fn, row.names = F, col.names = T)
