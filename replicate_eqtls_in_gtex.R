### this script collects eqtl status in different gtex tissues

library(argparser)
source('io_util.R')
source('arg_util.R')
library(corrplot)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtls with tissue-wise stats", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.per_tissue_stats.txt")
args <- add_argument(args, "-tis", help="tissues separated by comma", default="Whole_Blood,Muscle_Skeletal,Thyroid,Skin_Sun_Exposed_Lower_leg,Testis")
args <- add_argument(args, "-lab", help="tissues labels by comma", default="Whole Blood,Skeletal Muscle,Thyroid,Skin - Sun Exposed,Testis")
args <- add_argument(args, "-o", help="output file", default="results/eqtl_replication_gtex.pdf")

argv = parse_args(args)
eqtl_fn = argv$eqtl
tissue_input = argv$tis
label_input = argv$lab
out_fn = argv$o

fdr_th = 0.05

tissues = parse_delimitted_param(tissue_input, delim = ',', rm.empty = T)
labels = parse_delimitted_param(label_input, delim = ',', rm.empty = T)

### read data
eqtl_df = read_df(eqtl_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")

### plot start
pdf(out_fn)

# replicate eqtls from a tissue to other tissues
crossmappable_replication_n = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
crossmappable_replication_total = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
noncrossmappable_replication_n = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
noncrossmappable_replication_total = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))

per_gene_crossmappable_replication_n = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
per_gene_crossmappable_replication_total = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
per_gene_noncrossmappable_replication_n = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))
per_gene_noncrossmappable_replication_total = matrix(NA, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))

for(tissue in tissues){
  # get significant eqts in the tissue
  tissue_fdr_col = sprintf('FDR_%s', tissue)
  tissue_eqtl_df = eqtl_df[!is.na(eqtl_df[,tissue_fdr_col]) & eqtl_df[,tissue_fdr_col] <= fdr_th, ]
  crossmapped_tissue_eqtl_df = tissue_eqtl_df[tissue_eqtl_df$cross_mappable_genes != "", ]
  noncrossmapped_tissue_eqtl_df = tissue_eqtl_df[tissue_eqtl_df$cross_mappable_genes == "", ]
  
  ### replicate in other tissues
  other_tissues = setdiff(tissues, tissue)
  for(other in other_tissues){
    other_p_col = sprintf('pvalue_%s', other)
    crossmapped_p = crossmapped_tissue_eqtl_df[,other_p_col]
    crossmapped_p = crossmapped_p[!is.na(crossmapped_p)]  # excluding not-tested pairs
    noncrossmapped_p = noncrossmapped_tissue_eqtl_df[,other_p_col]
    noncrossmapped_p = noncrossmapped_p[!is.na(noncrossmapped_p)]
    
    ### plot p-value histograms in other tissues
    hist_breaks = seq(from = 0, to = 1, by = 0.05)
    crossmapped_col = rgb(0,0,1,1/4)
    noncrossmapped_col = rgb(1,0,0,1/4)
    h1 = hist(crossmapped_p, breaks = hist_breaks, plot = F)
    h2 = hist(noncrossmapped_p, breaks = hist_breaks, plot = F)
    h1$counts=h1$counts/sum(h1$counts)
    h2$counts=h2$counts/sum(h2$counts)
    plot( h1, col=crossmapped_col, 
          xlim=c(0,1), ylim = c(0,1),
          main = sprintf('pvalue : %s -> %s', tissue, other),
          xlab = 'p', ylab = 'Probability')
    plot( h2, col=noncrossmapped_col, xlim=c(0,1), add=T)
    legend('topright', 
           legend = c('Cross-mappable', 'Non-cross-mappable'), 
           fill = c(crossmapped_col, noncrossmapped_col), 
           bg = rgb(1,1,1,0.5))
    
    ### compute the number and fraction of replicated eqtls at fdr<=0.05
    other_fdr = p.adjust(c(crossmapped_p, noncrossmapped_p), method = 'BH')
    crossmappable_replication_total[tissue, other] = length(crossmapped_p)
    noncrossmappable_replication_total[tissue, other] = length(noncrossmapped_p)
    crossmappable_replication_n[tissue, other] = sum(other_fdr[1:length(crossmapped_p)] <= fdr_th)
    noncrossmappable_replication_n[tissue, other] = sum(other_fdr[length(crossmapped_p) + (1:length(noncrossmapped_p))] <= fdr_th)
    
    ### using best snp per gene, compute the number replicated eqtls at fdr<=0.05
    tissue_p_col = sprintf('pvalue_%s', tissue)
    other_p_col = sprintf('pvalue_%s', other)
    other_fdr_col = sprintf('FDR_%s', other)
    common_eqtl_df = tissue_eqtl_df[!is.na(tissue_eqtl_df[,other_p_col]), c('snps','gene',tissue_p_col, other_p_col, 'cross_mappable_genes')]
    common_eqtl_df = common_eqtl_df[order(common_eqtl_df[,tissue_p_col]), ]
    per_gene_common_eqtl_df = common_eqtl_df[! duplicated(common_eqtl_df$gene), ]
    per_gene_common_eqtl_df[,other_fdr_col] = p.adjust(per_gene_common_eqtl_df[,other_p_col], method='BH')
    crossmapped_per_gene_common_eqtl_df = per_gene_common_eqtl_df[per_gene_common_eqtl_df$cross_mappable_genes != "", ]
    noncrossmapped_per_gene_common_eqtl_df = per_gene_common_eqtl_df[per_gene_common_eqtl_df$cross_mappable_genes == "", ]
    per_gene_crossmappable_replication_total[tissue, other] = nrow(crossmapped_per_gene_common_eqtl_df)
    per_gene_noncrossmappable_replication_total[tissue, other] = nrow(noncrossmapped_per_gene_common_eqtl_df)
    per_gene_crossmappable_replication_n[tissue, other] = sum(crossmapped_per_gene_common_eqtl_df[,other_fdr_col] <= fdr_th)
    per_gene_noncrossmappable_replication_n[tissue, other] = sum(noncrossmapped_per_gene_common_eqtl_df[,other_fdr_col] <= fdr_th)
    
    ### qq-plot
    crossmapped_col = rgb(0,0,1)
    noncrossmapped_col = rgb(1,0,0)
    get_qq_points <- function(x){
      y <- sort(x)
      unif_x = (1:length(y)) / length(y)
      return(list(x=unif_x, y=y))
    }
    
    crossmapped_qq = get_qq_points(crossmapped_per_gene_common_eqtl_df[,tissue_p_col])
    noncrossmapped_qq = get_qq_points(noncrossmapped_per_gene_common_eqtl_df[,other_p_col])
    
    xmax = -log10(min(crossmapped_qq$x, noncrossmapped_qq$x))
    ymax = -log10(min(crossmapped_qq$y, noncrossmapped_qq$y))
    
    plot( -log10(crossmapped_qq$x), -log10(crossmapped_qq$y),
          col=crossmapped_col, 
          pch = 19,
          type='l', lty=1, 
          xlim=c(0,xmax), ylim = c(0,ymax),
          main = sprintf('QQ(using best SNP per gene) : %s -> %s', tissue, other),
          xlab = expression(paste("-log"["10"]~"(uniform p)")), ylab = expression(paste("-log"["10"]~"(p)")))
    points( -log10(crossmapped_qq$x), -log10(crossmapped_qq$y), pch=19, col=crossmapped_col)
    lines( -log10(noncrossmapped_qq$x), -log10(noncrossmapped_qq$y), col=noncrossmapped_col)
    points( -log10(noncrossmapped_qq$x), -log10(noncrossmapped_qq$y), pch=19, col=noncrossmapped_col)
    abline(a=0, b=1, lty=2)
    legend('topleft',
           legend = c('Cross-mappable', 'Non-cross-mappable'), 
           col = c(crossmapped_col, noncrossmapped_col), 
           lty = 1,
           pch = 19,
           bg = rgb(1,1,1,0.5))
  }
}

### plot fraction of replicated eqtls in cross-mappable and non-cross-mappable group
crossmappable_replication_frac = crossmappable_replication_n / crossmappable_replication_total
noncrossmappable_replication_frac = noncrossmappable_replication_n / noncrossmappable_replication_total

diag(crossmappable_replication_frac) <- 1
diag(noncrossmappable_replication_frac) <- 1

dimnames(crossmappable_replication_frac) = list(labels, labels)
dimnames(noncrossmappable_replication_frac) = list(labels, labels)

col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
corrplot(crossmappable_replication_frac, 
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         cl.lim = c(0, 1),
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction: cross-mappable"))

corrplot(noncrossmappable_replication_frac, 
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         cl.lim = c(0, 1),
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction: non-cross-mappable"))

corrplot(crossmappable_replication_frac - noncrossmappable_replication_frac,
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction: cross-mappable - non-cross-mappable"))


### using best snp per gene, plot fraction of replicated eqtls in cross-mappable and non-cross-mappable group
per_gene_crossmappable_replication_frac = per_gene_crossmappable_replication_n / per_gene_crossmappable_replication_total
per_gene_noncrossmappable_replication_frac = per_gene_noncrossmappable_replication_n / per_gene_noncrossmappable_replication_total

diag(per_gene_crossmappable_replication_frac) <- 1
diag(per_gene_noncrossmappable_replication_frac) <- 1

dimnames(per_gene_crossmappable_replication_frac) = list(labels, labels)
dimnames(per_gene_noncrossmappable_replication_frac) = list(labels, labels)

col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
corrplot(per_gene_crossmappable_replication_frac, 
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         cl.lim = c(0, 1),
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction (using best SNP per gene): cross-mappable"))

corrplot(per_gene_noncrossmappable_replication_frac, 
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         cl.lim = c(0, 1),
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction (using best SNP per gene): non-cross-mappable"))

corrplot(per_gene_crossmappable_replication_frac - per_gene_noncrossmappable_replication_frac,
         col = col(200),
         is.corr=T, 
         order = 'alphabet',
         method = "color", 
         addCoef.col = "black",
         number.cex = .8,
         number.digits = 2,
         diag = F,
         addgrid.col = "black",
         tl.col = 'black',
         title = sprintf("replicated fraction (using best SNP per gene): cross-mappable - non-cross-mappable"))


### plot end
dev.off()

### save data
rdata_fn = sprintf('%s.RData', out_fn)
save(crossmappable_replication_n, crossmappable_replication_total, crossmappable_replication_frac,
     noncrossmappable_replication_n, noncrossmappable_replication_total, noncrossmappable_replication_frac,
     per_gene_crossmappable_replication_n, per_gene_crossmappable_replication_total, per_gene_crossmappable_replication_frac,
     per_gene_noncrossmappable_replication_n, per_gene_noncrossmappable_replication_total, per_gene_noncrossmappable_replication_frac,
     tissues, labels, file = rdata_fn)
