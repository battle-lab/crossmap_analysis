# violin plot of correlation - crossmap vs non-crossmap

suppressMessages(library(argparser))
suppressMessages(source('io_util.R'))
suppressMessages(library(vioplot))

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-expr", help="expression data file", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/corrected_expression/Thyroid.v7.corrected.txt")
args <- add_argument(args, "-cross", help="cross-mappability file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt")
args <- add_argument(args, "-mincross", help="minimum cross-mappability", default = 100)
args <- add_argument(args, "-annot", help="gene annotation file (.txt)", default = '/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, "-map", help="mappability file", default = '/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/avg_mappability_Exon_UTR.txt')
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/crossmap_corr_violin.pdf")

argv = parse_args(args)
expr_fn = argv$expr
crossmap_fn = argv$cross
min_crossmap = argv$mincross
gene_annot_fn = argv$annot
mappability_fn = argv$map
plt_fn = argv$o

### read data
expr_df = read_df(expr_fn, sep = '\t', header = T, row.names = T)
expr_genes = rownames(expr_df)

gene_annot_df = read_df(gene_annot_fn, sep = '\t', header = T, row.names = T)
mappability_df = read_df(mappability_fn, sep = '\t', header = F, row.names = T)
crossmap_df = read_df(crossmap_fn, sep = '\t', header = F, row.names = F)

### take genes with non-NA mappability only
expr_genes_mappability = mappability_df[expr_genes, 1]
expr_genes = expr_genes[is.finite(expr_genes_mappability)]
expr_df = expr_df[expr_genes, ]

filter_crossmap_by_genes <- function(crossmap_df, incl.genes){
  # crossmap_df : cross-mappability data frame where first two columns are genes
  gene_fac = as.factor(crossmap_df[,1])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene1_included = gene_levels_included[as.integer(gene_fac)]
  
  gene_fac = as.factor(crossmap_df[,2])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene2_included = gene_levels_included[as.integer(gene_fac)]
  
  return(crossmap_df[gene1_included & gene2_included, , drop=F])
}

filter_crossmap_by_crossmappability <- function(crossmap_df, min.crossmap=-Inf, max.crossmap=Inf){
  # 3rd column contains cross-mappability
  if (min.crossmap > -Inf)
    crossmap_df = crossmap_df[crossmap_df[,3]>= min.crossmap, ]
  if (max.crossmap < Inf)
    crossmap_df = crossmap_df[crossmap_df[,3] <= max.crossmap, ]
  
  return(crossmap_df)  
}

get_crossmap_index <- function(cross_df, genes, threshold=-Inf){
  stopifnot(length(setdiff(c(cross_df[,1], cross_df[,2]), genes)) == 0)
  if(threshold > -Inf)
    cross_df = filter_crossmap_by_crossmappability(cross_df, min.crossmap = min_crossmap)
  
  genes1 = factor(cross_df[,1], levels = genes)
  genes2 = factor(cross_df[,2], levels = genes)
  
  crossmap_min_idx = pmin(as.integer(genes1), as.integer(genes2))
  crossmap_max_idx = pmax(as.integer(genes1), as.integer(genes2))
  crossmap_idx_in_matrix = unique((crossmap_min_idx-1) * length(genes) + crossmap_max_idx)
  
  return(crossmap_idx_in_matrix)  # returns integer index vector
}

get_non_crossmap_index <- function(cross_df, genes){
  lt = lower.tri(diag(length(genes)), diag = F)
  crossmap_idx_in_matrix = get_crossmap_index(cross_df, genes)
  lt[crossmap_idx_in_matrix] = F
  return(lt)  # returns logical index matrix
}

### compute pairwise correlation and split by cross-mappability
expr.cor = abs(cor(t(expr_df)))
tissue_crossmap_df = filter_crossmap_by_genes(crossmap_df, expr_genes)

crossmap_idx_in_matrix = get_crossmap_index(tissue_crossmap_df, expr_genes)
r_cross = expr.cor[crossmap_idx_in_matrix]

thresholded_crossmap_idx_in_matrix = get_crossmap_index(tissue_crossmap_df, expr_genes, threshold = min_crossmap)
r_cross_thresholded = expr.cor[thresholded_crossmap_idx_in_matrix]

non_crossmap_idx_in_matrix = get_non_crossmap_index(tissue_crossmap_df, expr_genes)
r_non_cross = expr.cor[non_crossmap_idx_in_matrix]

### compute pairwise correlation and split by cross-mappability for protein coding genes
is_protein_coding_genes = gene_annot_df[expr_genes, 'gene_type'] == 'protein_coding'
protein_coding_genes = expr_genes[is_protein_coding_genes]
coding.expr.cor = expr.cor[protein_coding_genes, protein_coding_genes]

tissue_crossmap_df = filter_crossmap_by_genes(tissue_crossmap_df, protein_coding_genes)

crossmap_idx_in_matrix = get_crossmap_index(tissue_crossmap_df, protein_coding_genes)
coding_r_cross = coding.expr.cor[crossmap_idx_in_matrix]

thresholded_crossmap_idx_in_matrix = get_crossmap_index(tissue_crossmap_df, protein_coding_genes, threshold = min_crossmap)
coding_r_cross_thresholded = coding.expr.cor[thresholded_crossmap_idx_in_matrix]

non_crossmap_idx_in_matrix = get_non_crossmap_index(tissue_crossmap_df, protein_coding_genes)
coding_r_non_cross = coding.expr.cor[non_crossmap_idx_in_matrix]

### wilcoxon test
#w_r = wilcox.test(x = r_cross, y = r_non_cross, alternative = 'g')
#w_r_th = wilcox.test(x = r_cross_thresholded, y = r_non_cross, alternative = 'g')

probs = c(0.5, 0.75, 0.95, 1)
quantlies_cross = quantile(r_cross, probs=probs)
quantlies_cross_thresholded = quantile(r_cross_thresholded, probs=probs)
quantlies_non_cross = quantile(r_non_cross, probs=probs)

quantlies_cross_str = paste(probs, format(quantlies_cross, scientific = F, digits = 3), sep=':', collapse = ', ')
quantlies_cross_thresholded_str= paste(probs, format(quantlies_cross_thresholded, scientific = F, digits = 3), sep=':', collapse = ', ')
quantlies_non_cross_str = paste(probs, format(quantlies_non_cross, scientific = F, digits = 3), sep=':', collapse = ', ')

coding_quantlies_cross = quantile(coding_r_cross, probs=probs)
coding_quantlies_cross_thresholded = quantile(coding_r_cross_thresholded, probs=probs)
coding_quantlies_non_cross = quantile(coding_r_non_cross, probs=probs)

coding_quantlies_cross_str = paste(probs, format(coding_quantlies_cross, scientific = F, digits = 3), sep=':', collapse = ', ')
coding_quantlies_cross_thresholded_str= paste(probs, format(coding_quantlies_cross_thresholded, scientific = F, digits = 3), sep=':', collapse = ', ')
coding_quantlies_non_cross_str = paste(probs, format(coding_quantlies_non_cross, scientific = F, digits = 3), sep=':', collapse = ', ')

### save data
data_fn = gsub(pattern = '.pdf$', replacement = '.pdf.RData', x = plt_fn)
save(r_cross, r_non_cross, r_cross_thresholded, coding_r_cross, coding_r_non_cross, coding_r_cross_thresholded, file = data_fn)

### plot
pdf(plt_fn)

# all genes
vioplot(r_cross, r_non_cross, names=c("Cross-mappable", "Non-cross-mappable"), col="gold")
title(main=paste0('Cross-quantiles: ', quantlies_cross_str, '\nNon-cross-quantiles: ', quantlies_non_cross_str), ylab = 'Correlation (|r|)')
vioplot(r_cross_thresholded, r_non_cross, names=c(paste0("Cross-mappable (>= ",min_crossmap," )"), "Non-cross-mappable"), col="gold")
title(main=paste0('Cross-quantiles (>=', min_crossmap, ') - ', quantlies_cross_thresholded_str, '\nNon-cross-quantiles: ', quantlies_non_cross_str), ylab = 'Correlation (|r|)')

# protein-coding genes
vioplot(coding_r_cross, coding_r_non_cross, names=c("Cross-mappable", "Non-cross-mappable"), col="gold")
title(main=paste0('Coding-cross-quantiles: ', coding_quantlies_cross_str, '\nCoding-non-cross-quantiles: ', coding_quantlies_non_cross_str), ylab = 'Correlation (|r|)')
vioplot(coding_r_cross_thresholded, coding_r_non_cross, names=c(paste0("Cross-mappable (>= ",min_crossmap,")"), "Non-cross-mappable"), col="gold")
title(main=paste0('Coding-cross-quantiles (>=', min_crossmap, '): ', coding_quantlies_cross_thresholded_str, '\nNon-cross-quantiles: ', coding_quantlies_non_cross_str), ylab = 'Correlation (|r|)')

dev.off()


