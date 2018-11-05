# filter eqtls based on cross-mappability withing 1Mb of the cross-mapping gene
# we get all snp-gene pairs with p<=1e-5
# but, to compute FDR, we need total number of tests after filtering.

library(argparser)
library(data.table)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file, all eqtls with p<=1e-5", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/trans_eqtl_cross_chr/Muscle_Skeletal_cross_chr_trans_eqtl_all_p_1e-5.txt")
args <- add_argument(args, "-gencode", help="gencode annotation file (txt format)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-snp_pos", help="file with snp positions", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process/combined_maf_0.05_repeat_masked.012.pos")
args <- add_argument(args, "-expr", help="expression data file, required to get genes tested", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Muscle_Skeletal.v7.normalized_expression.txt")
args <- add_argument(args, "-crossmap", help="crossmapping file", default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt')
args <- add_argument(args, "-d", help="distance threshold for cross-mappability", default=1e6)
args <- add_argument(args, "-o", help="annotated eqtl file", default="results/eqtl_crossmap_filtered.txt")

argv = parse_args(args)
eqtl_fn = argv$eqtl
gencode_fn = argv$gencode
snp_pos_fn = argv$snp_pos
expr_fn = argv$expr
crossmap_fn = argv$crossmap
d_th = argv$d
out_fn = argv$o

### read data
eqtl_df = read_df(eqtl_fn, sep = '\t', header = T, row.names = F)
gencode_df = read_df(gencode_fn, sep = '\t', header = T, row.names = F)
snp_pos_df = read_df(snp_pos_fn, sep = '\t', header = F, row.names = F)
expr_df = read_df(expr_fn, sep = '\t', header = T, row.names = T)
crossmap_df = read_df(crossmap_fn, sep = '\t', header = F, row.names = F)

### fix snp chr --- should have 'chr' prefix
if( ! startsWith(as.character(snp_pos_df[1,1]), 'chr'))
  snp_pos_df[,1] = paste0('chr', snp_pos_df[,1])

### which genes are tested
tested_genes = rownames(expr_df)

### function: fast acess genes location (crh, tss)
gene_loc_env = new.env(hash = T)
tss_values =  as.integer(apply(gencode_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))
gencode_df$tss = tss_values
tmp <- mapply(function(g, chr, tss){gene_loc_env[[g]] <<- list(chr=chr, tss=tss); return()}, gencode_df$gene_id, gencode_df$chr, gencode_df$tss)
get_gene_location <- function(g){
  return(list(chr=gene_loc_env[[g]]$chr, tss=gene_loc_env[[g]]$tss))
}

### function: given chr, get annotation of all genes in the given chr
chromosomes = unique(gencode_df$chr)
gene_anot_by_chr = lapply(chromosomes, function(chr) {gencode_df[gencode_df$chr==chr,,drop=F]})
names(gene_anot_by_chr) = chromosomes
get_gene_annot_in_chr <- function(chr){
  return(gene_anot_by_chr[[chr]])
}

### function: given a position, get nearby genes
### note: the genes do not have to be in the list of test genes
get_nearby_genes <- function(chr, pos, d=1e6){
  chr_annot_df = get_gene_annot_in_chr(chr)
  left_pos = pos - d
  right_pos = pos + d
  near_annot_df = chr_annot_df[chr_annot_df$tss >= left_pos & chr_annot_df$tss <=right_pos, ]
  return(near_annot_df$gene_id)
}

### for each gene, get cross-mappable genes from the given gene in different chr
### note: the cross-mappable genes must be in the list of test genes
crossmap_per_gene = tapply(crossmap_df$V2, crossmap_df$V1, c)
trans_crossmap_per_gene_env = new.env(hash = T)
tmp = lapply(names(crossmap_per_gene), function(g){
  cg = intersect(crossmap_per_gene[[g]], tested_genes) # cross-map gene must be tested
  g_chr = get_gene_location(g)$chr
  cg_chr = sapply(cg, function(x) get_gene_location(x)$chr)
  trans_cg = cg[cg_chr != g_chr]
  trans_crossmap_per_gene_env[[g]] <<- trans_cg
  return()
})

get_trans_crossmap_genes <- function(g){
  return(trans_crossmap_per_gene_env[[g]])
}

### function: given a snp location, 
### find which cross-mappable trans genes were tests with the snp
get_cross_test_per_snp <- function(chr, pos, d=1e6){
  genes = get_nearby_genes(chr, pos, d)
  cross_genes_list <- lapply(genes, get_trans_crossmap_genes)
  cross_genes = unique(unlist(cross_genes_list))
  return(cross_genes)
}

get_n_cross_test_per_snp <- function(chr, pos, d=1e6){
  return(length(get_cross_test_per_snp(chr, pos, d)))
}

### function: compute the total number of trans test made before filtering
get_n_total_test_without_cross_filtering <- function(snp_pos_df, tested_genes){
  snp_count_per_chr =  tapply(snp_pos_df[,1],  snp_pos_df[,1], length)
  test_genes_chr = sapply(tested_genes, function(g) get_gene_location(g)$chr)
  gene_count_per_chr = tapply(test_genes_chr, test_genes_chr, length)
  test_count_per_chr = sum(gene_count_per_chr) - gene_count_per_chr
  test_count_per_chr = test_count_per_chr[names(snp_count_per_chr)]
  total_test_count = sum(as.numeric(snp_count_per_chr) * test_count_per_chr)
  return(total_test_count)
}

### function: filter eqtl results for cross-mapping
filter_crossmap_eqtls <- function(eqtl_df, d=1e6){
  get_snp_loc_from_id <- function(snpid){
    parts = strsplit(snpid, split = '_')[[1]]
    chr = ifelse(startsWith(parts[1], prefix = 'chr'), parts[1], paste0('chr', parts[1]))
    pos = as.integer(parts[2])
    return(list(chr = chr, pos = pos))
  }
  
  is_crossmappable = mapply(function(snp_id, gene){
    snp_loc = get_snp_loc_from_id(snp_id)
    cg = get_cross_test_per_snp(chr = snp_loc$chr, pos = snp_loc$pos, d = d_th)
    return(gene %in% cg)
  }, eqtl_df$snps, eqtl_df$gene)
  
  filtered_eqtl_df = eqtl_df[!is_crossmappable, ,drop=F]
  return(filtered_eqtl_df)
}

#### compute total number of cross-mappable test
n_cross_test_per_snp <- mapply(get_n_cross_test_per_snp, snp_pos_df[,1], snp_pos_df[,2], rep(d_th, nrow(snp_pos_df)), SIMPLIFY = F)
total_n_cross_test <- sum(unlist(n_cross_test_per_snp))
total_test_before_filter <- get_n_total_test_without_cross_filtering(snp_pos_df, tested_genes)
total_test_after_filter <- total_test_before_filter - total_n_cross_test

### filter eqtls and update FDR
crossmap_filtered_eqtl_df <- filter_crossmap_eqtls(eqtl_df, d = d_th)
crossmap_filtered_eqtl_df$FDR <- p.adjust(crossmap_filtered_eqtl_df$pvalue, method = 'BH', n = total_test_after_filter)

### write output
write_df(crossmap_filtered_eqtl_df, file = out_fn, sep = "\t", row.names = F, col.names = T)
n_test_out_fn = gsub(pattern = '.txt$', replacement = '_ntest.txt', x = out_fn)
write_df(total_test_after_filter, file = n_test_out_fn, row.names = F, col.names = F)
