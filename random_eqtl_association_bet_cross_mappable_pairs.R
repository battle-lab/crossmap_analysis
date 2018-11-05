suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(vioplot))
suppressMessages(library(parallel))

args <- arg_parser('program')
args <- add_argument(args, '-meqtl', 
                     help='matrix eqtl directory',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls/Whole_Blood')
args <- add_argument(args, '-meta', 
                     help='matrix eqtl metadata file',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/matrix_eqtls_metadata/meqtl_meta_Whole_Blood.txt')
args <- add_argument(args, '-cis', 
                     help='best cis snp per gene',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/best_cis_snp_per_gene_1e6/best_cis_snp_Whole_Blood.txt')
args <- add_argument(args, '-cross', 
                     help='cross-mappability file',
                     default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt')
args <- add_argument(args, '-expr',
                     help='gene annotation file (txt)',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt')
args <- add_argument(args, '-map', 
                     help='gene-mappability file',
                     default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/avg_mappability_Exon_UTR.txt')
args <- add_argument(args, '-annot', 
                     help='gene annotation file',
                     default='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, '-n',
                     help='number of random pairs',
                     default=20)
args <- add_argument(args, '-sym',
                     help='is cross-mappability file symmetric',
                     default='FALSE')
args <- add_argument(args, '-q',
                     help='number of pairs in the top quantile',
                     default=10e3)
args <- add_argument(args, '-core',
                     help='number of cores',
                     default=1)
args <- add_argument(args, '-o',
                     help='output file. related data saved in a file with suffix .RData.',
                     default='results/random_eqtl_association.pdf')


argv <- parse_args(args)
matrix_eqtl_results_dir = argv$meqtl
matrix_eqtl_metadata_fn = argv$meta
best_cis_snp_fn = argv$cis
cross_map_fn = argv$cross
expr_fn = argv$expr
gene_map_fn = argv$map
gene_annot_fn = argv$annot
N0 = argv$n
is_symmetric_crossmap = as.logical(argv$sym)
n_pair_top_quantile = argv$q
n_cores = argv$core
plt_fn <- argv$o

N <- round(N0 * 1.1)  # take 10% extra eqtls, as some genes may not have any snp nearby

if(is.na(is_symmetric_crossmap) || is.null(is_symmetric_crossmap))
  is_symmetric_crossmap = FALSE

p_min = 1e-50
set.seed(101) # to reproduce results

### rad expression data
expr_df = read.table(expr_fn, sep ='\t', header = T, quote = "", comment.char = "", stringsAsFactors = F, row.names = 1)
dim(expr_df)
expr_df[1:5,1:5]

### read gene mappability file
gene_map_df = read.table(gene_map_fn, header = F, sep = '\t', row.names = 1, quote = "", comment.char = "")
dim(gene_map_df)
head(gene_map_df)

### read gene annotation file and compute tss
gene_annot_df = read.table(gene_annot_fn, header = T, sep = '\t', row.names = 1, quote = "", comment.char = "", stringsAsFactors = F)
dim(gene_annot_df)
head(gene_annot_df)

tss =  as.integer(apply(gene_annot_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))
gene_annot_df$tss = as.integer(tss)

### read matrix-eqtl metadata
meta_df = read.table(matrix_eqtl_metadata_fn, header = T, sep = '\t', quote = "", comment.char = "", stringsAsFactors = F)
meta_df$chr = as.character(meta_df$chr)
head(meta_df)

### read best cis snp per gene
best_cis_snp_df = read.table(best_cis_snp_fn, header = T, sep = '\t', quote = "", comment.char = "", stringsAsFactors = F)
rownames(best_cis_snp_df) = best_cis_snp_df$gene
head(best_cis_snp_df)

### read cross-map data
cross_map = fread(input = cross_map_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
dim(cross_map)
head(cross_map)

### keep only genes in the expression matrix
genes = rownames(expr_df)
col1_in_expr = cross_map$V1 %in% genes
col2_in_expr = cross_map$V2 %in% genes
cross_map = cross_map[col1_in_expr & col2_in_expr,]

### keep gene pairs from different chromosomes only
chr1 = gene_annot_df[cross_map$V1, 'chr']
chr2 = gene_annot_df[cross_map$V2, 'chr']
cross_map = cross_map[chr1!=chr2, ]

### get list of crossmapped-genes for each gene
if(!is_symmetric_crossmap){
  cross_map_per_gene = tapply(cross_map$V2, INDEX = cross_map$V1, FUN = c)
} else{
  cross_map_per_gene = tapply(c(cross_map$V2, cross_map$V1), INDEX = c(cross_map$V1, cross_map$V2), FUN = c)
}

### find genes in expression file with mappability NA
na_genes = intersect(rownames(gene_map_df)[is.na(gene_map_df$V2)], genes)
non_na_genes = intersect(rownames(gene_map_df)[!is.na(gene_map_df$V2)], genes)


######################### functions to compute cross-mapped eqtl association#########################
get_meqtl_files <- function(chr, pos, d){
  chr = gsub('chr', '', chr)
  chr_meta = meta_df[meta_df$chr == chr, ]
  left_pos = pos - d
  right_pos = pos + d
  outside = (chr_meta$start_pos<left_pos & chr_meta$end_pos<left_pos) | (chr_meta$start_pos>right_pos & chr_meta$end_pos>right_pos)
  chr_meta = chr_meta[!outside, ]
  return(chr_meta$file)
}

# file_meqtl_df_queue <-list()
# max_q_len = 5

get_meqtl_results <- function(fn){
  # # return if already in memory
  # if(fn %in% names(file_meqtl_df_queue))
  #   return(file_meqtl_df_queue[[fn]])
  # 
  # # if queue is already full, remove the oldest file
  # if(length(file_meqtl_df_queue) >= max_q_len){
  #   file_meqtl_df_queue <<-  file_meqtl_df_queue[-1]
  # }
  # 
  # load and insert the file
  stopifnot(file.exists(fn))
  load(fn)
  # ## add snp chr and snp pos
  # meqtl_snps_levels = levels(meqtl_df$snps)
  # meqtl_snps_levels_chr = sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) gsub(pattern = 'chr', replacement = '', x = parts[1], ignore.case = T) )
  # meqtl_df$snp_chr = meqtl_snps_levels_chr[as.integer(meqtl_df$snps)]
  # 
  # meqtl_snps_levels_pos = as.integer(sapply(strsplit(meqtl_snps_levels, split='_'), function(parts) parts[2]))
  # meqtl_df$snp_pos = meqtl_snps_levels_pos[as.integer(meqtl_df$snps)]
  # 
  # file_meqtl_df_queue[[fn]] <<- meqtl_df
  # return(file_meqtl_df_queue[[fn]])
  return(meqtl_df)
}

# get_best_snp <- function(g1){
#   # debug: g1="ENSG00000135297.11"
#   g1_chr_tss = gene_annot_df[g1, c('chr','tss')]
#   g1_chr = as.character(g1_chr_tss[1])
#   g1_tss = as.integer(g1_chr_tss[2])
#   g1_files <- get_meqtl_files(chr = g1_chr, pos = g1_tss, d = d_gene_snp)
#   
#   candidate_meqtl_results <- NULL
#   for(g1_fn in g1_files){
#     print(g1_fn)
#     res_meqtl_fn = sprintf('%s/%s', matrix_eqtl_results_dir, g1_fn)
#     res_meqtl_df = get_meqtl_results(fn = res_meqtl_fn)
#     candidates = res_meqtl_df[res_meqtl_df$gene == g1, ]
#     candidates = candidates[candidates$snp_pos>=(g1_tss-d_gene_snp) & candidates$snp_pos<=(g1_tss+d_gene_snp),  , drop=F]
#     if(nrow(candidates) > 0)
#       candidate_meqtl_results = rbind(candidate_meqtl_results, candidates)
#     rm(res_meqtl_df, candidates)
#   }
#   
#   if(is.null(candidate_meqtl_results))
#     return()
#   
#   best_snp_idx = which.min(candidate_meqtl_results$pvalue)
#   best_snp = as.character(candidate_meqtl_results$snps[best_snp_idx])
#   rm(candidate_meqtl_results)
#   gc(reset = T)
#   return(best_snp)
# }

get_best_snp <- function(g1){
  best_snp = best_cis_snp_df[g1, 'snps']
  if(is.na(best_snp))
    return()
  return(best_snp)
}

get_eqtl_association <- function(s, g){
  snp_chr = gsub(pattern = 'chr', replacement = '', x = strsplit(s, split='_')[[1]][1], ignore.case = T)
  snp_pos = as.integer(strsplit(s, split='_')[[1]][2])
  snp_file = get_meqtl_files(chr = snp_chr, pos = snp_pos, d = 0)
  stopifnot(length(snp_file)==1)
  meqtl_df = get_meqtl_results(sprintf('%s/%s', matrix_eqtl_results_dir, snp_file))
  association = meqtl_df[meqtl_df$snps == s & meqtl_df$gene == g, ]
  rm(meqtl_df)
  gc(reset = T)
  return(association)
}

get_crossmap_eqtl_association <- function(g1, g2){
  ## find the best snp associated with g1 (withing given distance)
  g1_best_snp = get_best_snp(g1 = g1)
  if(length(g1_best_snp)==0){
    return(c(beta=NA, p=NA))
  }
  
  ### return association between best snp and g2
  cross_assoc = get_eqtl_association(s = g1_best_snp, g = g2)
  return(list(beta=cross_assoc$beta, p=cross_assoc$pvalue))
}

### given a sampled crossmap dataframe, 
### return first N rows with at least one snp near gene 1
filter_crossmap_without_snp_near_gene1 <- function(sample_crossmap_df, maxrow=1000){
  best_snps = lapply(as.character(sample_crossmap_df[,1]), FUN = get_best_snp)
  has_best_snps = !sapply(best_snps, is.null)
  cum_count = cumsum(has_best_snps)
  to_take = has_best_snps
  to_take[cum_count > maxrow] = FALSE
  return(sample_crossmap_df[to_take, ])
}

######################### end of functions to compute crossmapped eqtl association #########################


### take random cross-map pairs
n_cross_map = nrow(cross_map)
sampled_cross_map = cross_map[sample(1:n_cross_map, size = N, replace = F, prob = cross_map$V3), ]
if(is_symmetric_crossmap){
  first_gene_idx = sample(1:2, size = nrow(sampled_cross_map), replace = T)
  gene1s = sampled_cross_map[,1]
  gene1s[which(first_gene_idx==2)] = sampled_cross_map[first_gene_idx==2,2]
  gene2s = sampled_cross_map[,2]
  gene2s[which(first_gene_idx==2)] = sampled_cross_map[first_gene_idx==2,1]
  sampled_cross_map[,1] = gene1s
  sampled_cross_map[,2] = gene2s
}

sampled_cross_map = filter_crossmap_without_snp_near_gene1(sampled_cross_map, maxrow = N0)

# ## sort sampled pairs to reduce i/o
# sampled_cross_map$g1_chr = gene_annot_df[sampled_cross_map[,1], 'chr', drop=F]
# sampled_cross_map$g1_tss = gene_annot_df[sampled_cross_map[,1], 'tss', drop=F]
# sampled_cross_map = sampled_cross_map[with(sampled_cross_map, order(g1_chr, g1_tss)), ]


### take random non-cross-mappable-pairs (avoid genes with NA mappability and pairs from same chr)
sampled_non_cross_map_genes1 = c()
sampled_non_cross_map_genes2 = c()
is_sampled = list()

while(length(sampled_non_cross_map_genes1) < N){
  g1 = sample(non_na_genes, size = 1, replace = F)
  non_cross_map_genes = setdiff(non_na_genes, g1)
  if(g1 %in% names(cross_map_per_gene))
    non_cross_map_genes = setdiff(non_cross_map_genes, cross_map_per_gene[[g1]])
  g2 = sample(non_cross_map_genes, size = 1)
  
  ### g1 and g2 has to be on different chr
  if(gene_annot_df[g1, 'chr'] == gene_annot_df[g2, 'chr'])
    next
  
  sampled_tag = paste0(g1,'|',g2)
  if(is_symmetric_crossmap && g1 > g2)
    sampled_tag = paste0(g2,'|',g1)
    
  if(!is.null(is_sampled[[sampled_tag]]))  # already sampled
    next
  
  sampled_non_cross_map_genes1 = c(sampled_non_cross_map_genes1, g1)
  sampled_non_cross_map_genes2 = c(sampled_non_cross_map_genes2, g2)
  is_sampled[[sampled_tag]] = TRUE
}

sampled_non_cross_map = data.frame(V1 = sampled_non_cross_map_genes1, V2 = sampled_non_cross_map_genes2, stringsAsFactors = F)

sampled_non_cross_map = filter_crossmap_without_snp_near_gene1(sampled_non_cross_map, maxrow = N0)


# ## sort sampled pairs to reduce i/o
# sampled_non_cross_map$g1_chr = gene_annot_df[sampled_non_cross_map[,1], 'chr', drop=F]
# sampled_non_cross_map$g1_tss = gene_annot_df[sampled_non_cross_map[,1], 'tss', drop=F]
# sampled_non_cross_map = sampled_non_cross_map[with(sampled_non_cross_map, order(g1_chr, g1_tss)), ]


### compute eqtl-association between sampled cross and non-cross mappable pairs
t1=Sys.time()
cross_map_assoc = mcmapply(FUN = get_crossmap_eqtl_association, 
                         sampled_cross_map$V1, 
                         sampled_cross_map$V2, 
                         mc.cores = n_cores, 
                         mc.cleanup = T)
gc(reset = T)
t2=Sys.time()
t2-t1

t1=Sys.time()
non_cross_map_assoc = mcmapply(FUN = get_crossmap_eqtl_association,
                             sampled_non_cross_map$V1,
                             sampled_non_cross_map$V2, 
                             mc.cores = n_cores, 
                             mc.cleanup = T)
gc(reset = T)
t2=Sys.time()
t2-t1


beta_cross = abs(as.numeric(cross_map_assoc["beta",]))
beta_cross = beta_cross[!is.na(beta_cross)]
p_cross = abs(as.numeric(cross_map_assoc["p",]))
p_cross = p_cross[!is.na(p_cross)]
p_cross_log = -log10(p_cross + p_min)

beta_non_cross = abs(as.numeric(non_cross_map_assoc["beta",]))
beta_non_cross = beta_non_cross[!is.na(beta_non_cross)]
p_non_cross = abs(as.numeric(non_cross_map_assoc["p",]))
p_non_cross = p_non_cross[!is.na(p_non_cross)]
p_non_cross_log = -log10(p_non_cross + p_min)


### test if correlation in cross mappable pairs is higher than non-cross-mappable pairs
w_beta = wilcox.test(x = beta_cross, y = beta_non_cross, alternative = 'g')
w_p = wilcox.test(x = p_cross, y = p_non_cross, alternative = 'l')

### store sampled correlations
crossmap_vs_noncrossmap_data = list()
crossmap_vs_noncrossmap_data[['sampled_cross_map']] = sampled_cross_map
crossmap_vs_noncrossmap_data[['sampled_non_cross_map']] = sampled_non_cross_map
crossmap_vs_noncrossmap_data[['beta_cross']] = beta_cross
crossmap_vs_noncrossmap_data[['p_cross']] = p_cross
crossmap_vs_noncrossmap_data[['beta_non_cross']] = beta_non_cross
crossmap_vs_noncrossmap_data[['p_non_cross']] = p_non_cross
crossmap_vs_noncrossmap_data[['w_beta']] = w_beta
crossmap_vs_noncrossmap_data[['w_p']] = w_p

save(crossmap_vs_noncrossmap_data, file = paste0(plt_fn, '.cross_nocross.RData'))

### sample from top quantile
top_quantile = sort(c((n_cross_map-c(n_pair_top_quantile))/n_cross_map))
top_qval = quantile(cross_map$V3, probs = top_quantile)
candidate_indexes = which(cross_map$V3 >= top_qval)
n_pairs_top_quantile <<- length(candidate_indexes)
sampled_cross_map_top_quantile = cross_map[sample(candidate_indexes, size = min(N,length(candidate_indexes)), replace = F), ]

if(is_symmetric_crossmap){
  first_gene_idx = sample(1:2, size = nrow(sampled_cross_map_top_quantile), replace = T)
  gene1s = sampled_cross_map_top_quantile[,1]
  gene1s[which(first_gene_idx==2)] = sampled_cross_map_top_quantile[first_gene_idx==2,2]
  gene2s = sampled_cross_map_top_quantile[,2]
  gene2s[which(first_gene_idx==2)] = sampled_cross_map_top_quantile[first_gene_idx==2,1]
  sampled_cross_map_top_quantile[,1] = gene1s
  sampled_cross_map_top_quantile[,2] = gene2s
}

sampled_cross_map_top_quantile = filter_crossmap_without_snp_near_gene1(sampled_cross_map_top_quantile, maxrow = N0)

### association in top quantiles
t1 = Sys.time()
top_quantile_cross_map_assoc = mcmapply(get_crossmap_eqtl_association, 
                                        sampled_cross_map_top_quantile$V1, 
                                        sampled_cross_map_top_quantile$V2,
                                        mc.cores = n_cores, 
                                        mc.cleanup = T)
gc(reset = T)
t2 = Sys.time()
t2 - t1

beta_top_quantile_cross = abs(as.numeric(top_quantile_cross_map_assoc["beta",]))
beta_top_quantile_cross = beta_top_quantile_cross[!is.na(beta_top_quantile_cross)]

p_top_quantile_cross = abs(as.numeric(top_quantile_cross_map_assoc["p",]))
p_top_quantile_cross = p_top_quantile_cross[!is.na(p_top_quantile_cross)]
p_top_quantile_cross_log = -log10(p_top_quantile_cross + p_min)

w_beta_top_quantile = wilcox.test(x = beta_top_quantile_cross, y = beta_non_cross, alternative = 'g')
w_p_top_quantile = wilcox.test(x = p_top_quantile_cross, y = p_non_cross, alternative = 'l')


### store top quantile's random association data
top_quantile_data = list()
top_quantile_data[['top_quantile_cross_map_assoc']] = top_quantile_cross_map_assoc
top_quantile_data[['top_quantile']] = top_quantile;
top_quantile_data[['top_qval']] = top_qval;
top_quantile_data[['n_pairs_top_quantile']] = n_pairs_top_quantile;
top_quantile_data[['w_beta_top_quantile']] = w_beta_top_quantile;
top_quantile_data[['w_p_top_quantile']] = w_p_top_quantile;


### save all data into file
save(crossmap_vs_noncrossmap_data, top_quantile_data, file = paste0(plt_fn, '.RData'))


### plot
pdf(plt_fn)
vioplot(beta_cross, beta_non_cross, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(ylab = 'beta')
vioplot(p_cross, p_non_cross, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p$p.value, scientific = T, digits = 2)), ylab = 'p')
vioplot(p_cross_log, p_non_cross_log, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p$p.value, scientific = T, digits = 2)), ylab = '-log10(p)')

vioplot(beta_top_quantile_cross, beta_non_cross, names=c("cross-mappable (top)", "non-cross-mappable"), col="gold")
title(ylab = 'beta')
vioplot(p_top_quantile_cross, p_non_cross, names=c("cross-mappable (top)", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross_top < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p_top_quantile$p.value, scientific = T, digits = 2)), ylab = 'p')
vioplot(p_top_quantile_cross_log, p_non_cross_log, names=c("cross-mappable (top)", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross_top < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p_top_quantile$p.value, scientific = T, digits = 2)), ylab = '-log10(p)')

dev.off()
