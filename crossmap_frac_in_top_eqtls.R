library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/Whole_Blood_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt")
args <- add_argument(args, "-gencode", help="gencode gene annotation file (txt format)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-rate", help="bacground crossmap rate", default=0.179109496062856)
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/crossmap_frac_in_top_eqtls.pdf")

argv = parse_args(args)
eqtl_fn = argv$eqtl
gencode_fn = argv$gencode
plt_fn = argv$o
bg_eqtl_crossmap_rate = argv$rate

### read eqtl file
eqtl_df = read.table(eqtl_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")
eqtl_df = eqtl_df[with(eqtl_df, order(FDR)),]

######### compute fraction of crossmap in top eQTLs #############
n_max_eqtl = nrow(eqtl_df)
if(n_max_eqtl < 4){
  first_ns = 1:n_max_eqtl
} else {
  first_ns = sort(unique( round(c(2^seq(2, log2(n_max_eqtl), by = 0.5), n_max_eqtl))))
}

crossmap_fractions = sapply(first_ns, function(first_n){
  frac = sum(eqtl_df$cross_mappable_genes[1:first_n] != "" & !is.na(eqtl_df$cross_mappable_genes[1:first_n]))/first_n
  return(frac)
})

data_fn = gsub(pattern = '.pdf$', replacement = '.pdf.RData', x = plt_fn)
save(first_ns, crossmap_fractions, bg_eqtl_crossmap_rate, file = data_fn)

### initialize pdf
pdf(plt_fn)

x = log2(first_ns)
y = crossmap_fractions
plot(x = x,
     y = y,
     col = 'black',
     lty = 1,
     type = 'l',
     xlim = c(0, max(x)),
     ylim = c(0, 1),
     xlab = expression(paste("log"["2"]~"(Number of top eQTLs)")),
     ylab = 'Proportion of cross-mappable eQTL pairs',
     main = 'Cross-mappable eQTL proportions in top eQTLs w/ cross-mappability')

points(x = x,
       y = y,
       col = 'black', 
       pch = 19)

if(bg_eqtl_crossmap_rate > 0)
  abline(h = bg_eqtl_crossmap_rate, col='red', lty=2)

dev.off()

