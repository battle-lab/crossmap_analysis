library(argparser)
library(MatrixEQTL)
library(parallel)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-files", help="shell script file patterns for matrix-eQTL results", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_crossmap/trans_eqtl/matrix_eqtls/Whole_Blood/meqtl_Whole_Blood*.RData")
args <- add_argument(args, "-filter", help="filter type: trans/map_trans/map_cross_trans", default="trans")
args <- add_argument(args, "-annot", help="gene annotation file (txt)", default='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, "-map", help="gene mappability file", default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/avg_mappability_Exon_UTR.txt')
args <- add_argument(args, "-crossmap", help="crossmapping file", default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength_symmetric_mean_sorted.txt')
args <- add_argument(args, "-thread", help="number of threads", default=1)
args <- add_argument(args, "-min_map", help="gene mappability threshold", default=0.8)
args <- add_argument(args, "-crossmap_d", help="threshold distance between snp and cross-mappable gene", default=500e3)
args <- add_argument(args, "-fdr", help="fdr thresholds, separated by comma", default='0.05,0.1,0.2')
args <- add_argument(args, "-o", help="out prefix", default="results/eqtl_filter")

argv = parse_args(args)
meqtl_files_pattern = argv$files
filter_type = argv$filter
gene_annot_fn = argv$annot
mappability_fn = argv$map
crossmap_fn = argv$crossmap
n_threads = argv$thread
gene_mappability_threshold = argv$min_map
crossmap_d = argv$crossmap_d
fdr_thresholds_input=argv$fdr
out_prefix = argv$o

# filter type must be one of the followings
stopifnot(filter_type %in% c('trans', 'map_trans', 'map_cross_trans'))

### get file names from pattern
meqtl_files = system(sprintf('ls %s', meqtl_files_pattern), intern = T)
stopifnot(is.null(attr(meqtl_files, 'status')))
stopifnot(length(meqtl_files)>0)

### parse fdr thresholds
parts = strsplit(fdr_thresholds_input, split = ',')[[1]]
parts = parts[sapply(parts, nchar)>0]
fdr_thresholds = as.numeric(parts)

### read annotation data
gene_annot_df = read_df(gene_annot_fn)
gene_annot_df$chr = gsub('chr', '', gene_annot_df$chr, ignore.case = T)

gene_map_df = NULL
if(filter_type %in% c('map_trans', 'map_cross_trans')){
  gene_map_df = read_df(mappability_fn, header = F)
  #colnames(gene_map_df) = 'mappability'
}
  
crossmap_df = NULL
if(filter_type %in% c('map_cross_trans')){
  crossmap_df = read_df(crossmap_fn, header = F, row.names = F)
  cross_map_per_gene = tapply(c(crossmap_df[,1], crossmap_df[,2]), INDEX = c(crossmap_df[,2], crossmap_df[,1]), FUN = unique)
  
  tss =  as.integer(apply(gene_annot_df, MARGIN = 1, FUN = function(row){
    ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
  }))
  names(tss) = rownames(gene_annot_df)
}

#### filters
filter_for_trans_associations <- function(meqtl_df){
  stopifnot(length(setdiff(levels(meqtl_df$gene), rownames(gene_annot_df)))==0 )
  meqtl_genes_levels = levels(meqtl_df$gene)
  meqtl_genes_levels_chr = gene_annot_df[meqtl_genes_levels, 'chr']
  meqtl_genes_chr = meqtl_genes_levels_chr[as.integer(meqtl_df$gene)]

  meqtl_snps_levels = levels(meqtl_df$snps)
  meqtl_snps_levels_chr = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) gsub(pattern = 'chr', replacement = '', x = parts[1], ignore.case = T) )
  meqtl_snps_chr = meqtl_snps_levels_chr[as.integer(meqtl_df$snps)]
  
  passed = meqtl_genes_chr!=meqtl_snps_chr
  return(passed)
}

filter_for_mappable_trans_associations <- function(meqtl_df){
  stopifnot(length(setdiff(levels(meqtl_df$gene), rownames(gene_annot_df)))==0 )
  meqtl_genes_levels = levels(meqtl_df$gene)
  meqtl_genes_levels_chr = gene_annot_df[meqtl_genes_levels, 'chr']
  meqtl_genes_chr = meqtl_genes_levels_chr[as.integer(meqtl_df$gene)]
  
  meqtl_snps_levels = levels(meqtl_df$snps)
  meqtl_snps_levels_chr = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) gsub(pattern = 'chr', replacement = '', x = parts[1], ignore.case = T) )
  meqtl_snps_chr = meqtl_snps_levels_chr[as.integer(meqtl_df$snps)]
  
  passed_cross_chr = meqtl_genes_chr!=meqtl_snps_chr
  
  meqtl_genes_levels_mappability = gene_map_df[meqtl_genes_levels, 1]
  meqtl_genes_mappability = meqtl_genes_levels_mappability[as.integer(meqtl_df$gene)]
  
  passed_gene_mappability = (!is.na(meqtl_genes_mappability)) & (meqtl_genes_mappability>=gene_mappability_threshold)
  
  passed = passed_cross_chr & passed_gene_mappability
  return(passed)
}

filter_for_mappable_crossmap_trans_associations <- function(meqtl_df){
  stopifnot(length(setdiff(levels(meqtl_df$gene), rownames(gene_annot_df)))==0 )
  meqtl_genes_levels = levels(meqtl_df$gene)
  meqtl_genes_levels_chr = gene_annot_df[meqtl_genes_levels, 'chr']
  meqtl_genes_chr = meqtl_genes_levels_chr[as.integer(meqtl_df$gene)]
  
  meqtl_snps_levels = levels(meqtl_df$snps)
  meqtl_snps_levels_chr = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) gsub(pattern = 'chr', replacement = '', x = parts[1], ignore.case = T) )
  meqtl_snps_chr = meqtl_snps_levels_chr[as.integer(meqtl_df$snps)]
  
  passed_cross_chr = meqtl_genes_chr!=meqtl_snps_chr
  
  meqtl_genes_levels_mappability = gene_map_df[meqtl_genes_levels, 1]
  meqtl_genes_mappability = meqtl_genes_levels_mappability[as.integer(meqtl_df$gene)]
  
  passed_gene_mappability = (!is.na(meqtl_genes_mappability)) & (meqtl_genes_mappability>=gene_mappability_threshold)
  
  meqtl_snps_levels_pos = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) parts[2])
  meqtl_snps_pos = meqtl_snps_levels_pos[as.integer(meqtl_df$snps)]
  
  meqtl_df$snp_chr = as.character(meqtl_snps_chr)
  meqtl_df$snp_pos = as.integer(meqtl_snps_pos)

  ############# approach #############
  # assumption: all genes are present id meqtl_df, all snps are from same chr and within a small region
  # - compute snp_region
  # - take all genes in snp_chr withing crossmap_d distance from snp_region
  # - for each gene in the region
  #   - sidx = get indexes snps within crossmap_d distance from its tss
  #   - pass[sidx] = meqtl_df[sidx, gene] %in% crosmapp(g)
  
  # compute snp_region
  meqtl_snp_chr = unique(meqtl_df$snp_chr)
  stopifnot(length(meqtl_snp_chr) == 1)
  
  leftmost_pos = min(meqtl_df$snp_pos) - crossmap_d
  rightmost_pos = min(meqtl_df$snp_pos) + crossmap_d
  
  # take all genes in snp_chr withing crossmap_d distance from snp_region
  nearby_genes = meqtl_genes_levels[meqtl_genes_levels_chr == meqtl_snp_chr]
  nearby_genes_tss = tss[nearby_genes]
  nearby_genes = nearby_genes[nearby_genes_tss >= leftmost_pos & nearby_genes_tss <= rightmost_pos]
  
  # - for each gene in the region
  #   - sidx = get indexes snps within crossmap_d distance from its tss
  #   - pass[sidx] = meqtl_df[sidx, gene] %in% crosmapp(g)
  
  passed_crossmap = rep(TRUE, nrow(meqtl_df))
  for(ng in nearby_genes){
    ng_crossmapped = cross_map_per_gene[[ng]]
    if(is.null(ng_crossmapped))
      next
    
    ng_tss = tss[ng]
    ng_left = ng_tss - crossmap_d
    ng_right = ng_tss + crossmap_d
    sidx = meqtl_df$snp_pos>= ng_left & meqtl_df$snp_pos<=ng_right
    meqtl_genes_levels_in_crossmap = meqtl_genes_levels %in% cross_map_per_gene[[ng]]
    meqtl_genes_in_crossmap = meqtl_genes_levels_in_crossmap[as.integer(meqtl_df[sidx, 'gene'])]
    passed_crossmap[sidx] = passed_crossmap[sidx] & !meqtl_genes_in_crossmap
  }
  
  passed = passed_cross_chr & passed_gene_mappability & passed_crossmap
  return(passed)
}

filter_meqtl_data <- function(meqtl_df, filter_type){
  if(filter_type == 'trans')
    return(filter_for_trans_associations(meqtl_df))
  if(filter_type == 'map_trans')
    return(filter_for_mappable_trans_associations(meqtl_df))
  if(filter_type == 'map_cross_trans')
    return(filter_for_mappable_crossmap_trans_associations(meqtl_df))
  return(NULL)
}


### compute fdr after filtering -- taking top pvalues
meqtl_test_info <- mclapply(meqtl_files, mc.cores = n_threads, FUN = function(meqtl_fn){
  print(meqtl_fn)
  if(exists('meqtl_df')){ rm(meqtl_df); gc(reset = T)}
  load(meqtl_fn)
  filter_pass_status <- filter_meqtl_data(meqtl_df, filter_type)
  small_p <- meqtl_df$p <= 1e-5
  selected_in_bh <- which(filter_pass_status & small_p)
  to_ret = list(meqtl_fn = meqtl_fn,
                n_tests = sum(filter_pass_status),
                meqtl_test = meqtl_df[selected_in_bh,,drop=F])
  
  rm(meqtl_df)
  gc(reset = T)
  return(to_ret)
})

### save intermediate results
meqtl_info_data_fn = sprintf("%s_meqtl_test_info.RData", out_prefix)
save(meqtl_test_info, file = meqtl_info_data_fn)


### compute fdr using gathered test pvalues
all_n_tests = sum(as.numeric(unlist(lapply(meqtl_test_info, function(x)x[['n_tests']]))))  # integer overflows
meqtl_results_df = do.call("rbind", lapply(meqtl_test_info, function(x)x[['meqtl_test']]))
meqtl_results_df$FDR = p.adjust(meqtl_results_df$pvalue, method = 'BH', n = all_n_tests)
meqtl_results_df = meqtl_results_df[order(meqtl_results_df$FDR), ]

### save results into files
meqtl_out_fn = sprintf("%s_all_p_1e-5.txt", out_prefix)
ntest_out_fn = sprintf("%s_ntest.txt", out_prefix)
write_df(meqtl_results_df, file = meqtl_out_fn, row.names = F, col.names = T)
write_df(all_n_tests, file = ntest_out_fn, row.names = F, col.names = F)

### save results after fdr correction
for(fdr in fdr_thresholds){
  fdr_eqtl_out_fn = sprintf("%s_fdr_%s.txt", out_prefix, fdr)
  egene_out_fn = sprintf("%s_fdr_%s_egene.txt", out_prefix, fdr)
  esnp_out_fn = sprintf("%s_fdr_%s_esnp.txt", out_prefix, fdr)
  fdr_eqtl_df = meqtl_results_df[meqtl_results_df$FDR <= fdr,]
  write_df(fdr_eqtl_df, file = fdr_eqtl_out_fn, row.names = F, col.names = T)
  write_df(unique(fdr_eqtl_df$gene), file = egene_out_fn, row.names = F, col.names = F)
  write_df(unique(fdr_eqtl_df$snps), file = esnp_out_fn, row.names = F, col.names = F)
}
