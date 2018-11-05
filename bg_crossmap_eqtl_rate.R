### this script computes the background rate of cross-mappable eqtls

library(argparser)
library(data.table)
source('io_util.R')

args <- arg_parser('program')
args <- add_argument(args, '-expr', 
                     help='expression files - comma separated',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Muscle_Skeletal.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Skin_Sun_Exposed_Lower_leg.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Testis.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Thyroid.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt')
args <- add_argument(args, '-geno', 
                     help='genotype directory',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process')
args <- add_argument(args, '-genosfx', 
                     help='genotype file suffix. note: prefix is chr name (e.g., chr3)',
                     default='_maf_0.05_repeat_masked.012.txt')
args <- add_argument(args, "-cross", 
                     help="cross-mappability file", 
                     default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt")
args <- add_argument(args, '-label',
                     help='Tissue labels (comma-separated)',
                     default='Muscle - Skeletal,Skin - Sun Exposed,Testis,Thyroid,Whole Blood')
args <- add_argument(args, "-annot", 
                     help="gencode annotation file (txt format)", 
                     default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/background_eqtl_crossmap.txt')

argv <- parse_args(args)
expr_fn_input = argv$expr
genotype_dir = argv$geno
genotype_suffix = argv$genosfx
crossmap_fn = argv$cross
labels_input = argv$label
gene_annot_fn = argv$annot
out_fn <- argv$o

chromosomes = paste0('chr', 1:22)

### parse inputs
parts = strsplit(expr_fn_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
expr_files = parts[is_valid_parts]

parts = strsplit(labels_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
tissue_labels = parts[is_valid_parts]

stopifnot(length(expr_files) == length(tissue_labels))

### read cross-map
crossmap_df = fread(input = crossmap_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
# S-G is cross-mappable, if reads from another gene g near S map to G
# given G, we need to find g such that crossmap(g->G) > 0
crossmap_per_gene = tapply(crossmap_df$V1, INDEX = crossmap_df$V2, FUN = c)

### read gencode annotation and compute TSS
gencode_df = read.table(gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")
rownames(gencode_df) = gencode_df$gene_id

#tss with rownames as gene id, columns: chromosome_name, transcription_start_site
tss_values =  as.integer(apply(gencode_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))

tss = data.frame(chromosome_name = gencode_df[gencode_df$gene_id, 'chr'],
                 transcription_start_site = tss_values, 
                 gene_id = gencode_df$gene_id,
                 row.names = gencode_df$gene_id, 
                 stringsAsFactors = F)


### read genotypes
chr_snps <- lapply(chromosomes, function(chr){
  print(paste('reading genotypes - ', chr))
  genotype_fn = sprintf("%s/%s%s", genotype_dir, chr, genotype_suffix)
  genotype_df = read_df(genotype_fn)
  snps = rownames(genotype_df)
  snps_pos = as.integer(sapply(strsplit(snps, '_'), function(parts) parts[2]))
  snps_df = data.frame(snp=snps, chr=chr, pos=snps_pos)
  return(snps_df)
})
names(chr_snps) = chromosomes

n_snps_per_chr = sapply(chr_snps, function(snp_df) nrow(snp_df))
n_possible_tests_per_gene_in_chr = sum(n_snps_per_chr) - n_snps_per_chr
n_possible_tests_per_gene_in_chr['chrX'] = sum(n_snps_per_chr)
n_possible_tests_per_gene_in_chr['chrY'] = sum(n_snps_per_chr)

### function to access tss fast
tss_gene_to_idx <- new.env(hash = T)
tss_genes = rownames(tss)
tmp <- sapply(1:length(tss_genes), function(idx) tss_gene_to_idx[[tss_genes[idx]]] <<- idx )

get_tss_entries <- function(genes, cols=NULL){
  tss_indexes = sapply(genes, function(g) tss_gene_to_idx[[g]])
  if(is.null(cols))
    return(tss[tss_indexes,,drop=F])
  return(tss[tss_indexes,cols,drop=F])
}

### function to access gencode fast
gencode_gene_to_idx <- new.env(hash = T)
gencode_genes = rownames(gencode_df)
tmp <- sapply(1:length(gencode_genes), function(idx) gencode_gene_to_idx[[gencode_genes[idx]]] <<- idx )

get_gencode_entries <- function(genes, cols=NULL){
  gencode_indexes = sapply(genes, function(g) gencode_gene_to_idx[[g]])
  if(is.null(cols))
    return(gencode_df[gencode_indexes,])
  return(gencode_df[gencode_indexes,cols])
}


### get snps near a gene
t1 = Sys.time()
snps_near_gene <- new.env(hash = T)
tmp = lapply(rownames(gencode_df), function(cg){
  d = 1e6
  #cross_tss_info = tss[cg,,drop=F]
  cross_tss_info = get_tss_entries(cg)
  cross_chr = cross_tss_info[cg, 'chromosome_name']
  cross_pos = cross_tss_info[cg, 'transcription_start_site']
  is_near_tss = chr_snps[[cross_chr]]$pos>= (cross_pos-d) & chr_snps[[cross_chr]]$pos <= (cross_pos+d)
  cross_snps = as.character(chr_snps[[cross_chr]][is_near_tss, 'snp'])
  snps_near_gene[[cg]] <<- cross_snps
  return()
})
t2 = Sys.time()
t2 - t1


bg_crossmap_rate = sapply(expr_files, function(expr_fn){
  print(expr_fn)
  expr_df = read_df(expr_fn, sep = '\t', row.names = T)
  genes = rownames(expr_df)
  
  get_trans_crossmap_count <- function(g){
    # get cross-mappable genes in different chr
    #cross_genes = intersect(crossmap_per_gene[[g]], genes)
    cross_genes = crossmap_per_gene[[g]]
    g_chr = gencode_df[g,'chr']
    cross_genes_chr = gencode_df[cross_genes, 'chr']
    cross_genes = cross_genes[cross_genes_chr != g_chr & cross_genes_chr %in% chromosomes]
    
    # get cross-mappinng snps
    # cross_snps = lapply(cross_genes, function(cg){
    #   cross_tss_info = tss[cg,,drop=F]
    #   cross_chr = tss[cg, 'chromosome_name']
    #   cross_pos = tss[cg, 'transcription_start_site']
    #   is_near_tss = chr_snps[[cross_chr]]$pos>= (cross_pos-d) & chr_snps[[cross_chr]]$pos <= (cross_pos+d)
    #   cross_snps = as.character(chr_snps[[cross_chr]][is_near_tss, 'snp'])
    #   return(cross_snps)
    # })
    # cross_snps = unique(unlist(cross_snps))
    
    cross_snps = lapply(cross_genes, function(cg) snps_near_gene[[cg]])
    cross_snps = unique(unlist(cross_snps))
    
    
    # get #trans-eqtl tests
    n_cross_trans_test = length(cross_snps)
    n_total_trans_test = as.numeric(n_possible_tests_per_gene_in_chr[g_chr])
    return(c(n_cross_trans_test=n_cross_trans_test, n_total_trans_test=n_total_trans_test))
  }
  
  trans_test_counts = sapply(genes, get_trans_crossmap_count)
  total_cross_tests = sum(trans_test_counts['n_cross_trans_test',])
  total_trans_tests = sum(trans_test_counts['n_total_trans_test',])
  bg_cross_eqtl_rate = total_cross_tests / total_trans_tests
  
  tmp_df = data.frame(fn= basename(expr_fn), total_trans_tests=total_trans_tests, total_cross_tests=total_cross_tests, bg_cross_eqtl_rate=bg_cross_eqtl_rate)
  write.table(tmp_df, file = paste0(out_fn, '.tmp'), row.names = F, col.names = T, quote = F, append = T)

  return(c(total_trans_tests=total_trans_tests, total_cross_tests=total_cross_tests, bg_cross_eqtl_rate=bg_cross_eqtl_rate))
})

bg_crossmap_rate_df = as.data.frame(t(bg_crossmap_rate))
bg_crossmap_rate_df$tissue = tissue_labels
bg_crossmap_rate_df$file = as.character(sapply(expr_files, basename))

write_df(bg_crossmap_rate_df, file = out_fn, row.names = F, col.names = T)
