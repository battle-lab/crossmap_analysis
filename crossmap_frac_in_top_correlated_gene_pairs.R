# this script explores the fraction of top cross-mappable pairs in top correlated pairs

suppressMessages(library(argparser))
suppressMessages(source('io_util.R'))

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-expr", help="expression data file", default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/corrected_expression/Thyroid.v7.corrected.txt")
args <- add_argument(args, "-cross", help="cross-mappability file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt")
args <- add_argument(args, "-overlap", help="overlap file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/positional_overlap.txt")
args <- add_argument(args, "-mincross", help="minimum cross-mappability", default = 1e-6)
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/crossmap_in_top_corr_gene_pairs_wo_overlap.pdf")

argv = parse_args(args)
expr_fn = argv$expr
crossmap_fn = argv$cross
overlap_fn = argv$overlap
min_crossmap = argv$mincross
plt_fn = argv$o

### read data
expr_df = read_df(expr_fn, sep = '\t', header = T, row.names = T)
expr_genes = rownames(expr_df)

overlap_df = read_df(overlap_fn, sep = '\t', header = T, row.names = F)

crossmap_df = read_df(crossmap_fn, sep = '\t', header = F, row.names = F)

#cross_mappability = tapply(crossmap_df$V2, INDEX = crossmap_df$V1, FUN = c)

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

tissue_crossmap_df = filter_crossmap_by_genes(crossmap_df, expr_genes)
tissue_crossmap_df = filter_crossmap_by_crossmappability(tissue_crossmap_df, min.crossmap = min_crossmap)

### process overlap data
tissue_overlap_df = filter_crossmap_by_genes(overlap_df, expr_genes)
overlap_genes1 = factor(tissue_overlap_df[,1], levels = expr_genes)
overlap_genes2 = factor(tissue_overlap_df[,2], levels = expr_genes)
overlap_min_idx = pmin(as.integer(overlap_genes1), as.integer(overlap_genes2))
overlap_max_idx = pmax(as.integer(overlap_genes1), as.integer(overlap_genes2))
overlap_idx_in_matrix = unique((overlap_min_idx-1) * length(expr_genes) + overlap_max_idx)

### compute the fraction of cross-mappable pairs at top
genes1 = factor(tissue_crossmap_df[,1], levels = expr_genes)
genes2 = factor(tissue_crossmap_df[,2], levels = expr_genes)

crossmap_min_idx = pmin(as.integer(genes1), as.integer(genes2))
crossmap_max_idx = pmax(as.integer(genes1), as.integer(genes2))

crossmap_idx_in_matrix = unique((crossmap_min_idx-1) * length(expr_genes) + crossmap_max_idx)

# non_crossmap_idx_in_matrix = lower.tri(diag(length(expr_genes)), diag = F)
# non_crossmap_idx_in_matrix[crossmap_idx_in_matrix] = F

n_genes = length(expr_genes)
expr.cor = abs(cor(t(expr_df)))

lt = lower.tri(diag(n_genes), diag = F)
lt[overlap_idx_in_matrix] = F # remove overlapped pairs
expr.cor.lower.tri = expr.cor[lt]
n_pairs = length(expr.cor.lower.tri)
expr.cor.lower.tri.rank = rank(-expr.cor.lower.tri, ties.method='random', na.last = T)

expr.cor.rank.mat = matrix(NA, nrow = n_genes, ncol = n_genes)
expr.cor.rank.mat[lt] = expr.cor.lower.tri.rank

#expr.cor.rank = rank(expr.cor, decreasing = T)
#crossmap.expr.cor.rank = sort(expr.cor.rank[crossmap_idx_in_matrix])

crossmap.expr.cor.rank = expr.cor.rank.mat[crossmap_idx_in_matrix]

### fraction of cross-mappabile pairs in top eQTLs
first_ns = sort(unique(round(c(2^seq(2, log2(n_pairs), by = 0.5), n_pairs))))
crossmap_fractions = sapply(first_ns, function(first_n){
  frac = sum(crossmap.expr.cor.rank <= first_n, na.rm = T)/ first_n  # na.rm=T to avoid NA due to overlap removal
  return(frac)
})

### save data
data_fn = gsub( pattern = '.pdf$', replacement = '.RData', plt_fn)
save(first_ns, crossmap_fractions, crossmap.expr.cor.rank, file = data_fn)

### plot
pdf(plt_fn)
x = log2(first_ns[-length(first_ns)])
y = crossmap_fractions[-length(first_ns)]
plot(x = x,
     y = y,
     col = 'coral3',
     lty = 1,
     type = 'l',
     xlim = c(0, max(x)),
     ylim = c(0, 1),
     xlab = "log2(Number of top correlated gene pairs)",
     ylab = 'Proportion of cross-mappable gene pairs',
     main = 'Cross-mappable gene-pair proportions in top correlated gene pairs')

points(x = x,
       y = y,
       col = 'coral3')

abline(h=crossmap_fractions[length(first_ns)], col='red')

dev.off()

