suppressMessages(library(argparser))
source('arg_util.R')

# args <- arg_parser("program");
# args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/prog/misc/cross_mappability/results/crossmap_frac_in_top_eqtls.pdf.RData")
# args <- add_argument(args, "-label", help="labels, comma separated", default="Whole Blood")
# args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined_crossmap_in_top_eqtls.pdf")

# args <- arg_parser("program");
# args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_protein_coding/crossmap_in_top_eqtl_p_1e-5/Muscle_Skeletal_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_protein_coding/crossmap_in_top_eqtl_p_1e-5/Skin_Sun_Exposed_Lower_leg_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_protein_coding/crossmap_in_top_eqtl_p_1e-5/Testis_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_protein_coding/crossmap_in_top_eqtl_p_1e-5/Thyroid_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_protein_coding/crossmap_in_top_eqtl_p_1e-5/Whole_Blood_crossmap_in_top_eqtls.pdf.RData")
# args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Sun Exposed,Testis,Thyroid,Whole Blood")
# args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined_crossmap_in_top_eqtls_protein_coding.pdf")


args <- arg_parser("program");
args <- add_argument(args, "-data", help="data files, comma separated", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_mappability_0.8/crossmap_in_top_eqtl_p_1e-5/Muscle_Skeletal_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_mappability_0.8/crossmap_in_top_eqtl_p_1e-5/Skin_Sun_Exposed_Lower_leg_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_mappability_0.8/crossmap_in_top_eqtl_p_1e-5/Testis_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_mappability_0.8/crossmap_in_top_eqtl_p_1e-5/Thyroid_crossmap_in_top_eqtls.pdf.RData,/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr_mappability_0.8/crossmap_in_top_eqtl_p_1e-5/Whole_Blood_crossmap_in_top_eqtls.pdf.RData")
args <- add_argument(args, "-label", help="labels, comma separated", default="Muscle - Skeletal,Skin - Sun Exposed,Testis,Thyroid,Whole Blood")
args <- add_argument(args, "-pt", help="font size in point", default=12)
args <- add_argument(args, "-width", help="plot width", default=7)
args <- add_argument(args, "-height", help="plot height", default=7)
args <- add_argument(args, "-margin", help="plot margin (bottom,left,up,right)", default=NA)
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/combined_crossmap_in_top_eqtls_mappability_0.8.pdf")


argv = parse_args(args)
data_files_input = argv$data
labels_input = argv$label
plt_font_size = argv$pt
plt_width = argv$width
plt_height = argv$height
plt_margin_input = argv$margin
plt_fn = argv$o

color_opts = c('darkolivegreen2', 'orchid', 'lightsalmon1', 'deepskyblue', 'slategray4', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'palegreen', 'bisque4', 'mediumpurple2', 'palevioletred3', 'bisque3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'tan1', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray', 'magenta',
               'navy', 'aquamarine', 'darkgoldenrod4', 'seagreen3', 'yellow3')

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
  if('first_ns' %in% ls()) rm(first_ns)
  if('bg_eqtl_crossmap_rate' %in% ls()) rm(bg_eqtl_crossmap_rate)
  load(data_fn)
  
  return(list(first_ns = first_ns,
              crossmap_fractions = crossmap_fractions,
              bg_eqtl_crossmap_rate = bg_eqtl_crossmap_rate))
})

### plot
max_x = ceiling(log2(max(unlist(lapply(combined_fractions, function(x) x$first_ns)))))
max_y = max(unlist(lapply(combined_fractions, function(x) max(x$crossmap_fractions))))
lty = 1

pdf(plt_fn, width = plt_width, height = plt_height, pointsize = plt_font_size)
if(!is.na(plt_margin_input) && nchar(plt_margin_input) > 0){
  margin = as.numeric(parse_delimitted_param(plt_margin_input, delim = ',', rm.empty = T))
  par(mai = margin)
}

for(di in 1:length(data_files)){
  x = log2(combined_fractions[[di]]$first_ns)
  y = combined_fractions[[di]]$crossmap_fractions
  bg_crossmap_rate = combined_fractions[[di]]$bg_eqtl_crossmap_rate
  tissue = tissue_labels[di]
  if(di == 1){
    plot(x = x,
         y = y,
         col = color_opts[tissue],
         lty = lty,
         type = 'l',
         xlim = c(0, max_x),
         ylim = c(0, max_y),
         xlab = expression(paste("log"["2"]~"(Number of top eQTLs)")),
         ylab = 'Proportion of cross-mappable eQTL pairs',
         main = 'Cross-mappable eQTL proportions in top eQTLs w/ cross-mappability')
  } else {
    points(x = x, y = y, col = color_opts[tissue], lty = lty, type = 'l')
  }
  
  points(x = x, y = y, col = color_opts[tissue], pch = 19)
  abline(h = bg_crossmap_rate, col = color_opts[tissue], lty = 3)
}

legend('topright', legend = tissue_labels, lwd = 1, lty = lty, col = color_opts[tissue_labels], pch = 19,  bg = rgb(1,1,1,0.5))


dev.off()
