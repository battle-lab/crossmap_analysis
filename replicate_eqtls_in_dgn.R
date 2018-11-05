### this script replicates eqtls in gtex-whole-blood in dgn

library(argparser)
source('io_util.R')
source('arg_util.R')
library(corrplot)
library(MatrixEQTL)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="whole blood eqtls in gtex", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/Whole_Blood_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt")
args <- add_argument(args, "-expr", help="dgn expr data", default="/work-zfs/abattle4/lab_data/dgn/data_used_for_eqtl_study/trans_data.txt")
args <- add_argument(args, "-geno", help="dgn genotype data", default="/work-zfs/abattle4/ashis/progdata/dgn/genotype/final_completed_genotype.txt")
args <- add_argument(args, "-cov", help="dgn covariate data", default="/work-zfs/abattle4/lab_data/dgn/covariates/Biological_and_hidden_factors.txt")
args <- add_argument(args, "-snp_annot", help="dgn snp annotation data", default="/work-zfs/abattle4/ashis/progdata/dgn/genotype/snp_annot.txt")
args <- add_argument(args, "-gene_annot", help="hg19 gene annotation data (txt)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-o", help="output file", default="results/eqtl_replication_dgn.txt")

argv = parse_args(args)
eqtl_fn = argv$eqtl
dgn_expr_fn = argv$expr
dgn_geno_fn = argv$geno
dgn_cov_fn = argv$cov
dgn_snp_annot_fn = argv$snp_annot
gene_annot_fn = argv$gene_annot
out_fn = argv$o

fdr_th = 0.05

### read data
eqtl_df = read_df(eqtl_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")
dgn_expr_df = read_df(dgn_expr_fn, sep = '\t', header = T, row.names = T, stringsAsFactors = F, quote = "")
dgn_geno_df = read_df(dgn_geno_fn, sep = '\t', header = T, row.names = T, stringsAsFactors = F, quote = "")
dgn_cov_df = read_df(dgn_cov_fn, sep = '\t', header = T, row.names = T, stringsAsFactors = F, quote = "")
dgn_snp_annot_df = read_df(dgn_snp_annot_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")
gene_annot_df = read_df(gene_annot_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")

### idential sample-order in dgn data
dgn_samples = intersect(rownames(dgn_expr_df), colnames(dgn_geno_df))
dgn_samples = intersect(dgn_samples, rownames(dgn_cov_df))
dgn_expr_df = dgn_expr_df[dgn_samples, ]
dgn_geno_df = dgn_geno_df[,dgn_samples]
dgn_cov_df = dgn_cov_df[dgn_samples, ]

### keep only sex, age, and genotype pcs in cov
dgn_cov_df = dgn_cov_df[, c("Sex", "PPAGE_at_interview", "geno PC1", "geno PC2", "geno PC3")]

### take eqtl at fdr<=0.05
eqtl_df = eqtl_df[eqtl_df$FDR<=fdr_th, ]
eGenes = unique(eqtl_df$gene)
eSNPs = unique(eqtl_df$snps)

### map ensembl geneids/snpids in eqtls to gene symbols / rsids in dgn
mapped_gene_annot_df = gene_annot_df[gene_annot_df$gene_id %in% eGenes, ]
mapped_gene_annot_df = mapped_gene_annot_df[mapped_gene_annot_df$gene_name %in% colnames(dgn_expr_df), ]
common_eqtl_df = merge(eqtl_df, mapped_gene_annot_df[,c('gene_id', 'gene_name')], by.x = 'gene', by.y = 'gene_id')
common_eqtl_df = merge(common_eqtl_df, dgn_snp_annot_df, by.x = c('snps_chr', 'snps_pos'), by.y = c('chr', 'pos'))
common_dgn_genes = common_eqtl_df$gene_name[!duplicated(common_eqtl_df$gene_name)]
common_gene_ids = common_eqtl_df$gene[!duplicated(common_eqtl_df$gene_name)]
common_dgn_snps = common_eqtl_df$snp[!duplicated(common_eqtl_df$snp)]
common_snp_ids = common_eqtl_df$snps[!duplicated(common_eqtl_df$snp)]

### prepare data for matrix-eqtl
common_dgn_expr_df = t(dgn_expr_df[dgn_samples, common_dgn_genes])
rownames(common_dgn_expr_df) = common_gene_ids
common_dgn_geno_df = dgn_geno_df[common_dgn_snps, dgn_samples]
rownames(common_dgn_geno_df) = common_snp_ids
common_dgn_cov_df = t(dgn_cov_df[dgn_samples, ])

meqtl_expr_slice = SlicedData$new(as.matrix(common_dgn_expr_df))
meqtl_snp_slice = SlicedData$new(as.matrix(common_dgn_geno_df))
meqtl_cov_slice = SlicedData$new(as.matrix(common_dgn_cov_df))

### run matrix-eQTL
me = Matrix_eQTL_engine(snps = meqtl_snp_slice,
                        gene = meqtl_expr_slice,
                        cvrt = meqtl_cov_slice,
                        output_file_name = NULL,
                        pvOutputThreshold = 1,
                        useModel = modelLINEAR, 
                        verbose = FALSE,
                        pvalue.hist = FALSE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE)

me$all$eqtls$snps = as.character(me$all$eqtls$snps)  # convert factor to chracter
me$all$eqtls$gene = as.character(me$all$eqtls$gene)  # convert factor to chracter

### merge new stats with given eqtls
dgn_eqtl_df = me$all$eqtls
common_eqtl_df = merge(common_eqtl_df, dgn_eqtl_df, by.x=c('snps','gene'), by.y=c('snps','gene'), suffixes = c('_gtex','_dgn'))
common_eqtl_df$FDR_dgn = p.adjust(common_eqtl_df$pvalue_dgn, method = 'BH')

### save
write_df(common_eqtl_df, file = out_fn)

### plot start
plt_fn = sprintf('%s.pdf', out_fn)
pdf(plt_fn)

### analyze dgn-replication
p_col = 'pvalue_dgn'
crossmapped_p = common_eqtl_df[common_eqtl_df$cross_mappable_genes != "",p_col]
crossmapped_p = crossmapped_p[!is.na(crossmapped_p)]  # excluding not-tested pairs
noncrossmapped_p = common_eqtl_df[common_eqtl_df$cross_mappable_genes == "",p_col]
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
      main = sprintf('pvalue : GTEx -> DGN'),
      xlab = 'p', ylab = 'Probability')
plot( h2, col=noncrossmapped_col, xlim=c(0,1), add=T)
legend('topright', 
       legend = c('Cross-mappable', 'Non-cross-mappable'), 
       fill = c(crossmapped_col, noncrossmapped_col), 
       bg = rgb(1,1,1,0.5))

### compute the number and fraction of replicated eqtls at fdr<=0.05
dgn_fdr = p.adjust(c(crossmapped_p, noncrossmapped_p), method = 'BH')
crossmappable_replication_total = length(crossmapped_p)
noncrossmappable_replication_total = length(noncrossmapped_p)
crossmappable_replication_n = sum(dgn_fdr[1:length(crossmapped_p)] <= fdr_th)
noncrossmappable_replication_n = sum(dgn_fdr[length(crossmapped_p) + (1:length(noncrossmapped_p))] <= fdr_th)
crossmappable_replication_frac = crossmappable_replication_n / crossmappable_replication_total
noncrossmappable_replication_frac = noncrossmappable_replication_n / noncrossmappable_replication_total

### plot
barplot_data = c(crossmappable_replication_frac, noncrossmappable_replication_frac)
names(barplot_data) = c('Cross-mappable', 'Non-cross-mappable') 
barplot(barplot_data, main = 'Replicated fraction: GTEx -> DGN', ylab = 'Replicated fraction')

### save plot data
rdata_fn = sprintf('%s.pdf.RData', out_fn)
save(barplot_data, 
     crossmappable_replication_n, crossmappable_replication_total, crossmappable_replication_frac,
     noncrossmappable_replication_n, noncrossmappable_replication_total, noncrossmappable_replication_frac,
     common_eqtl_df, file = rdata_fn)

### plot end
dev.off()

############################################################################################
############### analysis using best snp per gene, to avoid different #snps due to LD #######
############################################################################################
load(rdata_fn)
plt2_fn = sprintf('%s.best_snp_per_gene.pdf', out_fn)
rdata2_fn = sprintf('%s.best_snp_per_gene.pdf.RData', out_fn)
  
fdr_th = 0.05
crossmapped_col = rgb(0,0,1)
noncrossmapped_col = rgb(1,0,0)


common_eqtl_df = common_eqtl_df[with(common_eqtl_df, order(FDR_gtex)), ]
per_gene_common_eqtl_df = common_eqtl_df[! duplicated(common_eqtl_df$gene), ]
per_gene_common_eqtl_df$FDR_dgn = p.adjust(per_gene_common_eqtl_df$pvalue_dgn, method='BH')
crossmapped_per_gene_common_eqtl_df = per_gene_common_eqtl_df[per_gene_common_eqtl_df$cross_mappable_genes != "", ]
noncrossmapped_per_gene_common_eqtl_df = per_gene_common_eqtl_df[per_gene_common_eqtl_df$cross_mappable_genes == "", ]
crossmapped_per_gene_replication_frac = sum(crossmapped_per_gene_common_eqtl_df$FDR_dgn <= fdr_th) / nrow(crossmapped_per_gene_common_eqtl_df)
noncrossmapped_per_gene_replication_frac = sum(noncrossmapped_per_gene_common_eqtl_df$FDR_dgn <= fdr_th) / nrow(noncrossmapped_per_gene_common_eqtl_df)

pdf(plt2_fn)

### replicated fraction comparison
barplot_data = c(crossmapped_per_gene_replication_frac, noncrossmapped_per_gene_replication_frac)
names(barplot_data) = c('Cross-mappable', 'Non-cross-mappable') 
barplot(barplot_data, main = 'Replicated fraction using best SNP per gene: GTEx -> DGN', ylab = 'Replicated fraction')


### qq-plot
get_qq_points <- function(x){
  y <- sort(x)
  unif_x = (1:length(y)) / length(y)
  return(list(x=unif_x, y=y))
}

crossmapped_qq = get_qq_points(crossmapped_per_gene_common_eqtl_df$pvalue_dgn)
noncrossmapped_qq = get_qq_points(noncrossmapped_per_gene_common_eqtl_df$pvalue_dgn)

p_lim = 1e-16
#xmax = -log10(p_lim)
xmax = ceiling(-log10(min(crossmapped_qq$x, noncrossmapped_qq$x)))
ymax = -log10(p_lim)
crossmapped_qq$y[crossmapped_qq$y<=p_lim] = p_lim
noncrossmapped_qq$y[noncrossmapped_qq$y<=p_lim] = p_lim

plot( -log10(crossmapped_qq$x), -log10(crossmapped_qq$y),
      col=crossmapped_col, 
      pch = 20,
      #type='l', 
      lty=1, 
      xlim=c(0,xmax), ylim = c(0,ymax),
      main = sprintf('QQ(pvalue) : GTEx -> DGN'),
      xlab = expression(paste("-log"["10"]~"(uniform p)")), ylab = expression(paste("-log"["10"]~"(p)")))
abline(a=0, b=1, lty=2)
points( -log10(crossmapped_qq$x), -log10(crossmapped_qq$y), pch=19, col=crossmapped_col)
#lines( -log10(noncrossmapped_qq$x), -log10(noncrossmapped_qq$y), col=noncrossmapped_col)
points( -log10(noncrossmapped_qq$x), -log10(noncrossmapped_qq$y), pch=19, col=noncrossmapped_col)
legend('topleft',
       legend = c('Cross-mappable', 'Non-cross-mappable'), 
       col = c(crossmapped_col, noncrossmapped_col), 
       lty = 1,
       pch = 19,
       bg = rgb(1,1,1,0.5))
dev.off()

save(common_eqtl_df, per_gene_common_eqtl_df, 
     crossmapped_per_gene_common_eqtl_df, noncrossmapped_per_gene_common_eqtl_df, 
     crossmapped_per_gene_replication_frac, noncrossmapped_per_gene_replication_frac, 
     barplot_data, crossmapped_qq, noncrossmapped_qq, 
     xmax, ymax, 
     crossmapped_col, noncrossmapped_col,
     file = rdata2_fn)
