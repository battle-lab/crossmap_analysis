suppressMessages(library(argparser))

args <- arg_parser("program");

# args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Muscle_Skeletal.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Skin_Sun_Exposed_Lower_leg.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Testis.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Thyroid.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Whole_Blood.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.RData")
# args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Lower Leg,Testis,Thyroid,Whole Blood")
# args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_100.pdf")

# args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Muscle_Skeletal.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Skin_Sun_Exposed_Lower_leg.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Testis.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Thyroid.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Whole_Blood.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.RData")
# args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Lower Leg,Testis,Thyroid,Whole Blood")
# args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_1e-6.pdf")

# args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Muscle_Skeletal.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Skin_Sun_Exposed_Lower_leg.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Testis.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Thyroid.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Whole_Blood.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.RData")
# args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Lower Leg,Testis,Thyroid,Whole Blood")
# args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_100.pdf")

args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Muscle_Skeletal.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Skin_Sun_Exposed_Lower_leg.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Testis.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Thyroid.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/crossmap_in_top_correlated_gene_pairs/Whole_Blood.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.RData")
args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Lower Leg,Testis,Thyroid,Whole Blood")
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined.v7.corrected_crossmap_in_top_correlated_gene_pairs_threshold_1e-6.pdf")


argv = parse_args(args)
data_files_input = argv$data
labels_input = argv$label
plt_fn = argv$o

color_opts = c('darkolivegreen2', 'orchid', 'lightsalmon1', 'deepskyblue', 'slategray4', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'palegreen', 'bisque4', 'mediumpurple2', 'palevioletred3', 'bisque3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'tan1', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray')

# parse input
parts = strsplit(data_files_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
data_files = parts[is_valid_parts]

parts = strsplit(labels_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
tissue_labels = parts[is_valid_parts]

stopifnot(length(data_files)==length(tissue_labels))

color_opts = color_opts[1:length(tissue_labels)]
names(color_opts) = tissue_labels


### compute fractions
combined_fractions = lapply( data_files, function(data_fn){
  if('crossmap_fractions' %in% ls()) rm(crossmap_fractions)
  if('crossmap.expr.cor.rank' %in% ls()) rm(crossmap.expr.cor.rank)
  if('first_ns' %in% ls()) rm(first_ns)
  load(data_fn)
  
  max_n = max(first_ns)
  first_ns = sort(unique( round(c(2^seq(2, log2(max_n), by = 0.5), max_n))))
  crossmap_fractions = sapply(first_ns, function(first_n){
    frac = sum(crossmap.expr.cor.rank <= first_n, na.rm = T)/ first_n
    return(frac)
  })
  return(list(first_ns = first_ns,
              crossmap_fractions = crossmap_fractions))
})

### plot
max_x = ceiling(log2(max(unlist(lapply(combined_fractions, function(x) x$first_ns)))))
max_y = max(unlist(lapply(combined_fractions, function(x) max(x$crossmap_fractions))))
lty = 1
pdf(plt_fn)
for(di in 1:length(data_files)){
  x = log2(combined_fractions[[di]]$first_ns)
  y = combined_fractions[[di]]$crossmap_fractions
  tissue = tissue_labels[di]
  if(di == 1){
    plot(x = x,
         y = y,
         col = color_opts[tissue],
         lty = lty,
         type = 'l',
         xlim = c(0, max_x),
         ylim = c(0, max_y),
         xlab = "log2(Number of top correlated gene pairs)",
         ylab = 'Proportion of cross-mappable gene pairs',
         main = 'Cross-mappable gene-pair proportions in top correlated gene pairs')
  } else {
    points(x = x, y = y, col = color_opts[tissue], lty = lty, type = 'l')
  }
  
  points(x = x, y = y, col = color_opts[tissue], pch = 19)
}

legend('topright', legend = tissue_labels, lwd = 1, lty = lty, col = color_opts[tissue_labels], pch = 19,  bg = rgb(1,1,1,0.5))


dev.off()
