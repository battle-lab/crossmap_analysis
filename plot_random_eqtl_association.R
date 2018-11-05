library(argparser)

args <- arg_parser('program')
args <- add_argument(args, '-assoc', 
                     help='random eQTL association data (*.RData)',
                     #default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc/Whole_Blood_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc/Muscle_Skeletal_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc/Thyroid_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc/Testis_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc/Skin_Sun_Exposed_Lower_leg_random_eqtl_1e6.pdf.RData')
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc_v2_20180425/Whole_Blood_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc_v2_20180425/Muscle_Skeletal_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc_v2_20180425/Thyroid_random_eqtl_1e6.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/random_eqtl_assoc_v2_20180425/Skin_Sun_Exposed_Lower_leg_random_eqtl_1e6.pdf.RData')
args <- add_argument(args, '-label',
                     help='Tissue labels (comma-separated)',
                     default='Whole Blood,Muscle - Skeletal,Thyroid,Skin - Sun Exposed')
args <- add_argument(args, '-o',
                     help='output file (pdf)',
                     default='results/random_eqtl_assoc.pdf')

argv <- parse_args(args)
assoc_data_fn_input = argv$assoc
labels_input = argv$label
plt_fn <- argv$o

fdr_threshold = seq(0.05, 0.5, 0.05)

color_opts = c('darkolivegreen2', 'orchid', 'lightsalmon1', 'deepskyblue', 'slategray4', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'palegreen', 'bisque4', 'mediumpurple2', 'palevioletred3', 'bisque3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'tan1', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray')

### parse inputs
parts = strsplit(assoc_data_fn_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
assoc_data_files = parts[is_valid_parts]

parts = strsplit(labels_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
tissue_labels = parts[is_valid_parts]

color_opts = color_opts[1:length(tissue_labels)]
names(color_opts) = tissue_labels

stopifnot(length(assoc_data_files) == length(tissue_labels))


### function to count number of eQTL hits
get_n_eqtl_hits <- function(p, fdr_threshold=seq(0.05, 0.5, 0.05)){
  fdr_vals = p.adjust(p, method='BH')
  counts = sapply(fdr_threshold, function(th) sum(fdr_vals<=th) )
  return(counts)
}

### process each file
n_eqtl_hits = lapply(assoc_data_files,  function(assoc_data_fn){
  ### load data
  load(assoc_data_fn)

  p_cross = crossmap_vs_noncrossmap_data$p_cross
  p_non_cross = crossmap_vs_noncrossmap_data$p_non_cross
  p_top_cross = as.numeric(top_quantile_data$top_quantile_cross_map_assoc["p",])
  
  # # for debugging
  # N = 960
  # p_cross = p_cross[!is.na(p_cross)][1:N]
  # p_non_cross = p_non_cross[!is.na(p_non_cross)][1:N]
  # p_top_cross = p_top_cross[!is.na(p_top_cross)][1:N]

  ### check if category class has same number of eqtls
  stopifnot(length(p_cross) == length(p_non_cross))
  stopifnot(length(p_cross) == length(p_top_cross))
  
  n_cross = get_n_eqtl_hits(p_cross, fdr_threshold = fdr_threshold)
  n_non_cross = get_n_eqtl_hits(p_non_cross, fdr_threshold = fdr_threshold)
  n_top_cross = get_n_eqtl_hits(p_top_cross, fdr_threshold = fdr_threshold)
  
  ### clear memory
  rm(crossmap_vs_noncrossmap_data, top_quantile_data)
  gc(reset = T)
  
  ### return 
  return(list(n_cross=n_cross, n_non_cross=n_non_cross, n_top_cross=n_top_cross))
})
names(n_eqtl_hits) = tissue_labels

### plot
max_hits = max(unlist(n_eqtl_hits))
y_lim = ceiling(max_hits/5) * 5
line_width = 1

pdf(plt_fn)

for(ti in 1:length(tissue_labels)){
  tissue = tissue_labels[ti]
  n_cross = n_eqtl_hits[[tissue]]$n_cross
  n_non_cross = n_eqtl_hits[[tissue]]$n_non_cross
  n_top_cross = n_eqtl_hits[[tissue]]$n_top_cross
  
  if(ti == 1){
    plot(x = fdr_threshold,
         y = n_non_cross,
         col = color_opts[tissue],
         lty = 3,
         type = 'l',
         lwd = line_width+1,
         ylim = c(0, y_lim),
         xlab = "FDR Threshold",
         ylab = 'Number of trans-eQTLs',
         main = '')
  } else {
    lines(x = fdr_threshold, y = n_non_cross, col = color_opts[tissue], lty = 3, type = 'l', lwd = line_width+1)
  }
  
  points(fdr_threshold, n_non_cross, col = color_opts[tissue], pch=0)
  
  lines(x = fdr_threshold, y = n_cross, col = color_opts[tissue], lty = 2, type = 'l', lwd = line_width)
  points(fdr_threshold, n_cross, col = color_opts[tissue], pch=19)
  
  lines(x = fdr_threshold, y = n_top_cross, col = color_opts[tissue], lty = 1, type = 'l', lwd = line_width)
  points(fdr_threshold, n_top_cross, col = color_opts[tissue], pch=17)
  
}

legends_labels = c('Cross-mappable (Top)', 'Cross-mappable', 'Non Cross-mappable', tissue_labels)
legends_lty = c(1, 2, 3, rep(1,length(tissue_labels)))
legends_pch = c(17, 19, 0, rep(20,length(tissue_labels)))
legends_col = c('black', 'black', 'black', color_opts[tissue_labels])
legend('topleft', legend = legends_labels, lwd = 1, lty = legends_lty, col = legends_col, pch = legends_pch,  bg = rgb(1,1,1,0.5))

dev.off()
