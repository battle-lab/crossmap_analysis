library(argparser)
library(data.table)
library(vioplot)

args <- arg_parser('program')
args <- add_argument(args, '-cross', 
                     help='cross-mappability file',
                     default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt')
args <- add_argument(args, '-expr',
                     help='gene annotation file (txt)',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt')
args <- add_argument(args, '-map', 
                     help='gene-mappability file',
                     default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/avg_mappability_Exon_UTR.txt')
args <- add_argument(args, '-n',
                     help='number of random pairs',
                     default=1000)
args <- add_argument(args, '-n0',
                     help='number of random pairs for cross-map vs non-cross-map comparison',
                     default=1000)
args <- add_argument(args, '-sym',
                     help='is cross-mappability file symmetric',
                     default='FALSE')
args <- add_argument(args, '-q',
                     help='quantile separator - comma separated (0 and 1 are added). empty string for default partition.',
                     default='')
args <- add_argument(args, '-sep',
                     help='cross mappability separator values - comma separated',
                     default='1,2,5,10,20,50,100,200,300,400,500,1e8')
args <- add_argument(args, '-o',
                     help='output file. related data saved in a file with suffix .RData.',
                     default='results/random_cor_crossmap.pdf')


argv <- parse_args(args)
cross_map_fn = argv$cross
expr_fn = argv$expr
gene_map_fn = argv$map
N = argv$n
N0 = argv$n0
is_symmetric_crossmap = as.logical(argv$sym)
quantile_separators_input = argv$q
crossmap_separators_input = argv$sep
plt_fn <- argv$o

stopifnot(N0>=N)  # #samples in each group <= #samples for cross-map vs non-cross-map 

if(is.na(is_symmetric_crossmap) || is.null(is_symmetric_crossmap))
  is_symmetric_crossmap = FALSE

p_min = 1e-50
set.seed(101) # to reproduce results

crossmap_separators_values = as.numeric(strsplit(crossmap_separators_input, split = ',')[[1]])
stopifnot(length(crossmap_separators_values)>0)

quantile_separators_values = as.numeric(strsplit(quantile_separators_input, split = ',')[[1]])
#stopifnot(length(quantile_separators_values)>0)

# expr_fn = "/scratch1/battle-fs1/ashis/progdata/gtex_v7/per_tissue/WholeBlood.txt"
# cross_map_fn = "/scratch0/battle-fs1/annotation/mappability/pairwise_conflict.txt"
# plt_fn = "/scratch1/battle-fs1/ashis/results/misc/cross_mappability/random_cor_WholeBlood.pdf"
# N = 1000   # number of random pairs

# expr_fn = "/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/WholeBlood.txt"
# cross_map_fn = "/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_cross_mappability_strength.txt"
# plt_fn = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/random_cor_WholeBlood.pdf"
# N = 1000   # number of random pairs
# min_avg_tpm = 1

### rad expression data
expr_df = read.table(expr_fn, sep ='\t', header = T, quote = "", comment.char = "", stringsAsFactors = F, row.names = 1)
dim(expr_df)
expr_df[1:5,1:5]

# expr_df = expr_df[,-1]
# expr_df[1:5,1:5]

# ### keep genes with minimum average expression
# mean_tpm = rowMeans(expr_df)
# expr_df = expr_df[mean_tpm>=min_avg_tpm, ]

### read gene mappability file
gene_map_df = read.table(gene_map_fn, header = F, sep = '\t', row.names = 1, quote = "", comment.char = "")
dim(gene_map_df)
head(gene_map_df)

### read cross-map data
cross_map = fread(input = cross_map_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
dim(cross_map)
head(cross_map)

### keep only genes in the expression matrix
genes = rownames(expr_df)
col1_in_expr = cross_map$V1 %in% genes
col2_in_expr = cross_map$V2 %in% genes
cross_map = cross_map[col1_in_expr & col2_in_expr,]

### get list of crossmapped-genes for each gene
if(!is_symmetric_crossmap){
  cross_map_per_gene = tapply(cross_map$V2, INDEX = cross_map$V1, FUN = c)
} else{
  cross_map_per_gene = tapply(c(cross_map$V2, cross_map$V1), INDEX = c(cross_map$V1, cross_map$V2), FUN = c)
}

### find genes in expression file with mappability NA
na_genes = intersect(rownames(gene_map_df)[is.na(gene_map_df$V2)], genes)
non_na_genes = intersect(rownames(gene_map_df)[!is.na(gene_map_df$V2)], genes)

### take random cross-map pairs
n_cross_map = nrow(cross_map)
sampled_cross_map = cross_map[sample(1:n_cross_map, size = N0, replace = F, prob = cross_map$V3), ]

# ### take random non-cross-mappable-pairs -- with one gene from sampled cross-mappable pairs
# cross_map_per_gene = tapply(cross_map$V2, INDEX = cross_map$V1, FUN = c)
# 
# sampled_non_cross_map_genes2 = sapply(sampled_cross_map$V1, function(g1){
#   non_cross_map_genes = setdiff(genes, c(g1, cross_map_per_gene[[g1]]))
#   sample(non_cross_map_genes, size = 1, replace = F)
# })
# 
# sampled_non_cross_map = data.frame(V1 = sampled_cross_map$V1, V2 = sampled_non_cross_map_genes2, stringsAsFactors = F)

### take random non-cross-mappable-pairs (avoid genes with NA mappability)
sampled_non_cross_map_genes1 = c()
sampled_non_cross_map_genes2 = c()
is_sampled = list()

while(length(sampled_non_cross_map_genes1) < N0){
  g1 = sample(non_na_genes, size = 1, replace = F)
  non_cross_map_genes = setdiff(non_na_genes, g1)
  if(g1 %in% names(cross_map_per_gene))
    non_cross_map_genes = setdiff(non_cross_map_genes, cross_map_per_gene[[g1]])
  g2 = sample(non_cross_map_genes, size = 1)
  
  if(is_symmetric_crossmap){
    if(g1>g2){
      tmp = g1
      g1 = g2
      g2 = tmp
    }
  }

  sampled_tag = paste0(g1,'|',g2)
  if(!is.null(is_sampled[[sampled_tag]]))  # already sampled
    next
    
  sampled_non_cross_map_genes1 = c(sampled_non_cross_map_genes1, g1)
  sampled_non_cross_map_genes2 = c(sampled_non_cross_map_genes2, g2)
  is_sampled[[sampled_tag]] = TRUE
}

sampled_non_cross_map = data.frame(V1 = sampled_non_cross_map_genes1, V2 = sampled_non_cross_map_genes2, stringsAsFactors = F)


### compute correlation between sampled cross and non-cross mappable pairs
cross_map_cor = mapply(FUN = function(g1, g2){  
                               test = cor.test(as.numeric(expr_df[g1,]), as.numeric(expr_df[g2,]))
                               return(list(r=test$estimate, p=test$p.value))
                             },
                      sampled_cross_map$V1,
                      sampled_cross_map$V2)


non_cross_map_cor = mapply(FUN = function(g1, g2){  
                            test = cor.test(as.numeric(expr_df[g1,]), as.numeric(expr_df[g2,]))
                            return(list(r=test$estimate, p=test$p.value))
                          },
                          sampled_non_cross_map$V1,
                          sampled_non_cross_map$V2)

r_cross = abs(as.numeric(cross_map_cor["r",]))
r_cross = r_cross[!is.na(r_cross)]
p_cross = abs(as.numeric(cross_map_cor["p",]))
p_cross = p_cross[!is.na(p_cross)]
p_cross_log = -log10(p_cross + p_min)

r_non_cross = abs(as.numeric(non_cross_map_cor["r",]))
r_non_cross = r_non_cross[!is.na(r_non_cross)]
p_non_cross = abs(as.numeric(non_cross_map_cor["p",]))
p_non_cross = p_non_cross[!is.na(p_non_cross)]
p_non_cross_log = -log10(p_non_cross + p_min)


### test if correlation in cross mappable pairs is higher than non-cross-mappable pairs
w_r = wilcox.test(x = r_cross, y = r_non_cross, alternative = 'g')
w_p = wilcox.test(x = p_cross, y = p_non_cross, alternative = 'l')

### store sampled correlations
crossmap_vs_noncrossmap_data = list()
crossmap_vs_noncrossmap_data[['sampled_cross_map']] = sampled_cross_map
crossmap_vs_noncrossmap_data[['sampled_non_cross_map']] = sampled_non_cross_map
crossmap_vs_noncrossmap_data[['r_cross']] = r_cross
crossmap_vs_noncrossmap_data[['p_cross']] = p_cross
crossmap_vs_noncrossmap_data[['r_non_cross']] = r_non_cross
crossmap_vs_noncrossmap_data[['p_non_cross']] = p_non_cross
crossmap_vs_noncrossmap_data[['w_r']] = w_r
crossmap_vs_noncrossmap_data[['w_p']] = w_p


### function to compute pairwise correlation between genes
pairwise_cor_func = function(g1, g2){  
  test = cor.test(as.numeric(expr_df[g1,]), as.numeric(expr_df[g2,]))
  return(list(r=test$estimate, p=test$p.value))
}

### samples pairs from quanties of cross-mappability values
if(length(quantile_separators_values) > 0){
  quantiles = sort(unique(c(quantile_separators_values, 1)))
} else {
  #quantiles = sort(c(2000/n_cross_map, seq(0.1,0.9,0.1), 0.95, 1, (n_cross_map-c(2e3,4e3,6e3,8e3,10e3,15e3,25e3,50e3,100e3))/n_cross_map))
  right_quantiles = c((n_cross_map-c(2e3,4e3,6e3,8e3,10e3,15e3,25e3,50e3,100e3))/n_cross_map,1)
  left_quantiles = c(2000/n_cross_map, seq(0.1,0.9,0.1), 0.95)
  left_quantiles = left_quantiles[left_quantiles<min(right_quantiles)]
  quantiles = sort(unique(c(left_quantiles, right_quantiles)))
}

qvals = quantile(cross_map$V3, probs = quantiles)
qvals_selected = !duplicated(qvals)
qvals = qvals[qvals_selected]
quantiles = quantiles[qvals_selected]

n_pairs_quantiles <- NULL
quantile_sampled_cross_map_list <- lapply(0:length(qvals), function(qidx){
  if(qidx==0){
    # sample from non-cross-mappable pairs
    sampled_cross_map = sampled_non_cross_map[1:N,]
    return(sampled_cross_map)
  } else {
    q1 = ifelse(qidx == 1, 0, qvals[qidx-1])
    q2 = qvals[qidx]
    candidate_indexes = which(cross_map$V3 > q1 & cross_map$V3 <= q2)
    n_pairs_quantiles <<- c(n_pairs_quantiles, length(candidate_indexes))
    sampled_cross_map = cross_map[sample(candidate_indexes, size = min(N,length(candidate_indexes)), replace = F), ]
  }
})
names(quantile_sampled_cross_map_list) <- as.character(c(0, quantiles))

quantile_cross_map_cor_list <- lapply(quantile_sampled_cross_map_list, function(quantile_sampled_cross_map){
  quantile_cross_map_cor = mapply(pairwise_cor_func, quantile_sampled_cross_map$V1, quantile_sampled_cross_map$V2)
  r = abs(as.numeric(quantile_cross_map_cor["r",]))
  r = r[!is.na(r)]
  p = abs(as.numeric(quantile_cross_map_cor["p",]))
  p = p[!is.na(p)]
  return(list(r=r, p=p))
})

quantile_cross_map_wilcox_list = lapply(quantile_cross_map_cor_list, function(cor_list){
  base_corr_list = quantile_cross_map_cor_list[[1]]  # first one is non-crossmap samples
  w_r = wilcox.test(x = cor_list[['r']], y = base_corr_list[['r']], alternative = 'g')
  w_p = wilcox.test(x = cor_list[['p']], y = base_corr_list[['p']], alternative = 'l')
  return(list(w_r=w_r, w_p=w_p))
})


### store quantile-group's random correlation data
quantile_group_data = list()
quantile_group_data[['quantile_sampled_cross_map_list']] = quantile_sampled_cross_map_list
quantile_group_data[['quantile_cross_map_cor_list']] = quantile_cross_map_cor_list;
quantile_group_data[['quantiles']] = quantiles;
quantile_group_data[['qvals']] = qvals;
quantile_group_data[['n_pairs_quantiles']] = n_pairs_quantiles;
quantile_group_data[['quantile_cross_map_wilcox_list']] = quantile_cross_map_wilcox_list;

### samples pairs from selected cross-mappability range
qvals = crossmap_separators_values
n_pairs_range <- NULL
range_sampled_cross_map_list <- lapply(0:length(qvals), function(qidx){
  if(qidx==0){ 
    sampled_cross_map = sampled_non_cross_map[1:N,]
    return(sampled_cross_map)
  } else {
    q1 = ifelse(qidx == 1, 0, qvals[qidx-1])
    q2 = qvals[qidx]
    candidate_indexes = which(cross_map$V3 > q1 & cross_map$V3 <= q2)
    n_pairs_range <<- c(n_pairs_range, length(candidate_indexes))
    sampled_cross_map = cross_map[sample(candidate_indexes, size = min(N, length(candidate_indexes)), replace = F), ]
  }
})
names(range_sampled_cross_map_list) <- as.character(c(0, qvals))

range_cross_map_cor_list <- lapply(range_sampled_cross_map_list, function(quantile_sampled_cross_map){
  quantile_cross_map_cor = mapply(pairwise_cor_func, quantile_sampled_cross_map$V1, quantile_sampled_cross_map$V2)
  r = abs(as.numeric(quantile_cross_map_cor["r",]))
  r = r[!is.na(r)]
  p = abs(as.numeric(quantile_cross_map_cor["p",]))
  p = p[!is.na(p)]
  return(list(r=r, p=p))
})

range_cross_map_wilcox_list = lapply(range_cross_map_cor_list, function(cor_list){
  base_corr_list = range_cross_map_cor_list[[1]]  # first one is non-crossmap samples
  w_r = wilcox.test(x = cor_list[['r']], y = base_corr_list[['r']], alternative = 'g')
  w_p = wilcox.test(x = cor_list[['p']], y = base_corr_list[['p']], alternative = 'l')
  return(list(w_r=w_r, w_p=w_p))
})

### store quantile-group's random correlation data
range_group_data = list()
range_group_data[['range_sampled_cross_map_list']] = range_sampled_cross_map_list
range_group_data[['range_cross_map_cor_list']] = range_cross_map_cor_list;
range_group_data[['crossmap_vals']] = qvals;
range_group_data[['n_pairs_range']] = n_pairs_range;
range_group_data[['range_cross_map_wilcox_list']] = range_cross_map_wilcox_list;

### save all data into file
save(crossmap_vs_noncrossmap_data, quantile_group_data, range_group_data, file = paste0(plt_fn, '.RData'))


### plot
pdf(plt_fn)
vioplot(r_cross, r_non_cross, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(main=paste0('R_cross > R_non_cross?\nWilcoxon rank sum p <= ', format(w_r$p.value, scientific = T, digits = 2)), ylab = 'Correlation (r)')
vioplot(p_cross, p_non_cross, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p$p.value, scientific = T, digits = 2)), ylab = 'p')
vioplot(p_cross_log, p_non_cross_log, names=c("cross-mappable", "non-cross-mappable"), col="gold")
title(main=paste0('P_cross < P_non_cross?\nWilcoxon rank sum p <= ', format(w_p$p.value, scientific = T, digits = 2)), ylab = '-log10(p)')

par(las=2)
par(mar=c(8,5,2,1))

qcl = quantile_cross_map_cor_list
qs = names(quantile_cross_map_cor_list)
qs = sapply(qs, substr, 1, 6)
ws = quantile_cross_map_wilcox_list

arglist = lapply(qcl, function(x) x$r)
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_quantiles ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability quantiles',
      ylab = 'Correlation (r)')

arglist = lapply(qcl, function(x) x$p)
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_quantiles ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability quantiles',
      ylab = 'p-value')

arglist = lapply(qcl, function(x) -log10(x$p + p_min))
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_quantiles ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability quantiles',
      ylab = '-log10(p-value)')


wp = sapply(ws, function(x) -log10(x$w_r$p.value + p_min))
names(wp) = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_quantiles ,')'))
barplot(wp, main="Wilcoxon p-values\n(compared to non-cross-mapping group)",
        names.arg=names(wp),
        ylab='-log10(p)')


qcl = range_cross_map_cor_list
qs = names(range_cross_map_cor_list)
qs = sapply(qs, substr, 1, 6)
ws = range_cross_map_wilcox_list

arglist = lapply(qcl, function(x) x$r)
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_range ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability range',
      ylab = 'Correlation (r)')

arglist = lapply(qcl, function(x) x$p)
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_range ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability range',
      ylab = 'p-value')

arglist = lapply(qcl, function(x) -log10(x$p + p_min))
names(arglist) = 'x'
arglist[['names']] = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_range ,')'))
arglist[['col']]="gold"
do.call(vioplot, args = arglist)
title(main=paste0('Correlation between gene pairs'),
      xlab = 'Cross-mappability range',
      ylab = '-log10(p-value)')

wp = sapply(ws, function(x) -log10(x$w_r$p.value + p_min))
names(wp) = c("0", paste0(qs[1:(length(qs)-1)], '-', qs[2:length(qs)], '(', n_pairs_range ,')'))
barplot(wp, main="Wilcoxon p-values\n(compared to non-cross-mapping group)",
        names.arg=names(wp),
        ylab='-log10(p)')

dev.off()
