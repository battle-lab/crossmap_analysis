library(argparser)
library(scales)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", 
                     help="eqtl replication in rsem", 
                     default="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/replication_in_rsem/rsem_eqtl/Muscle_Skeletal_trans_eqtl_replication_in_rsem_fdr_0.05_allvar.txt")
args <- add_argument(args, "-lab", 
                     help="Dataset label used in plots", 
                     default='RSEM')
args <- add_argument(args, "-o", 
                     help="plot file", 
                     default="results/rsem_pvalue_categorized_by_crossmap.pdf")

argv = parse_args(args)
eqtl_fn = argv$eqtl
dataset_label = argv$lab
out_fn = argv$o


### read data
eqtl_df = read_df(eqtl_fn, row.names = F)
eqtl_df$col = c('red', 'black')[ as.numeric(eqtl_df$cross_mappable_genes=="") + 1]
eqtl_df$col2 = c('darkgreen', 'gray85')[ as.numeric(eqtl_df$gene_type=="pseudogene") + 1]
eqtl_df$col3 = c('red', 'darkgreen')[ as.numeric(eqtl_df$cross_mappable_genes=="") + 1]
eqtl_df$pch = c(4,1)[ as.numeric(eqtl_df$gene_type=="pseudogene") + 1]

pdf(out_fn)

### scatter plot gtex vs dgn pvalue
plot(-log10(eqtl_df$linreg_struct.pval), -log10(eqtl_df$pvalue_dgn),
     main = "",
     xlab = bquote("-log"["10"]~"(GTEx p)"),
     ylab = bquote("-log"["10"]~"("~.(dataset_label)~"p)"), 
     col = alpha(eqtl_df$col, 0.4))
abline(a=0, b=1, lty='dashed')
abline(h = 0, lty='dashed')
legend('topleft', pch = 1, legend = c('Cross-mappable', 'Not cross-mappable'), col =  c('red', 'black'), bg = alpha('white',0.5))

plot(-log10(eqtl_df$linreg_struct.pval), -log10(eqtl_df$pvalue_dgn),
     main = "",
     xlab = bquote("-log"["10"]~"(GTEx p)"),
     ylab = bquote("-log"["10"]~"("~.(dataset_label)~"p)"), 
     col = alpha(eqtl_df$col2, 0.4))
abline(a=0, b=1, lty='dashed')
abline(h = 0, lty='dashed')
legend('topleft', pch = 1, legend = c('Pseudogene target', 'Other target'), col =  c('gray85', 'darkgreen'), bg = alpha('white',0.5))

plot(-log10(eqtl_df$linreg_struct.pval), -log10(eqtl_df$pvalue_dgn),
     main = "",
     xlab = bquote("-log"["10"]~"(RNA-SeQC p)"),
     ylab = bquote("-log"["10"]~"("~.(dataset_label)~"p)"), 
     col = alpha(eqtl_df$col3, 0.4),
     pch = eqtl_df$pch)
abline(a=0, b=1, lty='dashed')
abline(h = 0, lty='dashed')
legend('topleft', pch = c(19,19,1,4), legend = c('Cross-mappable', 'Not cross-mappable', 'Pseudogene target', 'Other target' ), col =  c('red', 'darkgreen', 'black', 'black'), bg = alpha('white',0.5))


dev.off()

