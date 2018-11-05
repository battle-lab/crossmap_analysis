# plots random-correlation for different cross-mappabilities

library(argparser)
library(data.table)
source('plot_util.R')

args <- arg_parser('program')
args <- add_argument(args, '-data', 
                     help='data for random correlations',
                     default='/Users/ashissaha/github/misc/cross_mappability/results/gtex_v7/random_corr/Whole_Blood_random_correlation_uncorrected.pdf.RData,/Users/ashissaha/github/misc/cross_mappability/results/gtex_v7_kmer_count_normalized/random_corr/Whole_Blood_random_correlation_uncorrected.pdf.RData')
args <- add_argument(args, '-label',
                     help='labels for plots (comma separated)',
                     default='Regular,Normalized')
args <- add_argument(args, "-pt", help="font size in point", default=12)
args <- add_argument(args, "-width", help="plot width", default=20)
args <- add_argument(args, "-height", help="plot height", default=8)
args <- add_argument(args, '-o',
                     help='output file (pdf)',
                     default='results/random_cor_comparison.pdf')

argv <- parse_args(args)
data_files_input = argv$data
labels_input = argv$label
plt_font_size = argv$pt
plt_width = argv$width
plt_height = argv$height
plt_fn <- argv$o

data_files = as.character(strsplit(data_files_input, split = ',')[[1]])
stopifnot(length(data_files)>0)

data_labels = as.character(strsplit(labels_input, split = ',')[[1]])
stopifnot(length(data_labels)>0)
stopifnot(length(data_files)==length(data_labels))

#color_opts = c('deepskyblue', 'lemonchiffon', 'plum', 'palegreen', 'pink', 'gray', 'cyan', 'purple', 'blue')
color_opts = c('darkolivegreen2', 'orchid', 'lightsalmon1', 'deepskyblue', 'slategray4', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'palegreen', 'bisque4', 'mediumpurple2', 'palevioletred3', 'bisque3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'tan1', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray', 'magenta',
               'navy', 'aquamarine', 'darkgoldenrod4', 'seagreen3', 'yellow3')

stopifnot(length(data_files)<=length(color_opts))
color_opts = color_opts[1:length(data_files)]
names(color_opts) = data_labels[1:length(data_labels)]

### load data
all_data = lapply(data_files, function(fn){
  load(fn)
  return(list(crossmap_vs_noncrossmap_data=crossmap_vs_noncrossmap_data,
              quantile_group_data=quantile_group_data))
})
names(all_data)=data_labels

### plot
pdf(plt_fn, width = plt_width, height = plt_height, pointsize = plt_font_size)
par(las=2)
par(mar=c(10,5,2,1))


### plot crossmap_vs_noncrossmap
arg_list = list()
arg_labels = c()
arg_colors = c()

# put non-cross-mapping corr first
for(lbl in data_labels){
  arg_list[[length(arg_list)+1]] = all_data[[lbl]]$crossmap_vs_noncrossmap_data$r_non_cross
  arg_labels[length(arg_labels)+1] = sprintf('No cross-mapping')
  arg_colors[length(arg_colors)+1] = color_opts[lbl]
}

# put cross-mapping corr first
for(lbl in data_labels){
  arg_list[[length(arg_list)+1]] = all_data[[lbl]]$crossmap_vs_noncrossmap_data$r_cross
  arg_labels[length(arg_labels)+1] = sprintf('Cross-mapping')
  arg_colors[length(arg_colors)+1] = color_opts[lbl]
}

names(arg_list) = 'x'
arg_list[['names']] = arg_labels
arg_list[['col']] = arg_colors
arg_list['pchMed'] = '.:'

do.call(colored_vioplot, args = arg_list)
title(main=paste0('Correlation between gene pairs'),
      xlab = '',
      ylab = 'Correlation (r)')
legend("topleft", fill = color_opts,
       legend = data_labels, box.lty=0, bg = rgb(1,1,1,0.5), horiz = T)






################################################################
### plot crossmap_vs_noncrossmap where tissues stay together ###
################################################################
arg_list = list()
arg_labels = c()
arg_colors = c()

for(lbl in data_labels){
  arg_list[[length(arg_list)+1]] = all_data[[lbl]]$crossmap_vs_noncrossmap_data$r_cross
  arg_labels[length(arg_labels)+1] = sprintf('Cross-mappable')
  arg_colors[length(arg_colors)+1] = color_opts[lbl]
  
  arg_list[[length(arg_list)+1]] = all_data[[lbl]]$crossmap_vs_noncrossmap_data$r_non_cross
  arg_labels[length(arg_labels)+1] = sprintf('Not cross-mappable')
  arg_colors[length(arg_colors)+1] = color_opts[lbl]
}

names(arg_list) = 'x'
arg_list[['names']] = arg_labels
arg_list[['col']] = arg_colors
arg_list['pchMed'] = '.:'

do.call(colored_vioplot, args = arg_list)
title(main=paste0('Correlation between gene pairs'),
      xlab = '',
      ylab = 'Correlation (|r|)')
legend("topleft", fill = color_opts,
       legend = data_labels, box.lty=0, bg = rgb(1,1,1,0.5), horiz = F)

################################################




### plot wilcoxon p-values for crossmap_vs_noncrossmap
wp = sapply(data_labels, function(lbl){-log10(all_data[[lbl]]$crossmap_vs_noncrossmap_data$w_r$p.value + 1e-100)})
barplot(wp, main="Wilcoxon p-values (R_cross > R_non_cross)",
        names.arg=names(wp),
        ylab=expression(-log[10]*p), 
        col = color_opts)



### plot quantiles : v1
arg_list = list()
arg_labels = c()
arg_colors = c()

n_groups = length(all_data[[1]]$quantile_group_data$quantiles)+1
for(lbl in data_labels){
  for(gr in 1:n_groups){
    arg_list[[length(arg_list)+1]] = all_data[[lbl]]$quantile_group_data$quantile_cross_map_cor_list[[gr]]$r
    axis_lbl = ""
    if(gr==1){
      axis_lbl = sprintf("#%d: No cross-mapping", gr)
    } else{
      q1 = ifelse(gr==2, 0, all_data[[lbl]]$quantile_group_data$quantiles[gr-2])
      q2 = all_data[[lbl]]$quantile_group_data$quantiles[gr-1]
      q_pairs = all_data[[lbl]]$quantile_group_data$n_pairs_quantiles[gr-1]
      axis_lbl = sprintf("#%d: %.4g-%.4g(%d)",gr, q1, q2, q_pairs)
    }
    
    arg_labels[length(arg_labels)+1] = axis_lbl
    arg_colors[length(arg_colors)+1] = color_opts[lbl]
  }
}


names(arg_list) = 'x'
arg_list[['names']] = arg_labels
arg_list[['col']] = arg_colors
arg_list['pchMed'] = '.:'

do.call(colored_vioplot, args = arg_list)
title(main=paste0('Correlation between gene pairs'),
      xlab = '',
      ylab = 'Correlation (r)')
legend("topleft", fill = color_opts,
       legend = data_labels, box.lty=0, bg = rgb(1,1,1,0.5), horiz = T)



### plot quantiles : v2
arg_list = list()
arg_labels = c()
arg_colors = c()

n_groups = length(all_data[[1]]$quantile_group_data$quantiles)+1
for(gr in 1:n_groups){
  for(lbl in data_labels){
    arg_list[[length(arg_list)+1]] = all_data[[lbl]]$quantile_group_data$quantile_cross_map_cor_list[[gr]]$r
    axis_lbl = ""
    if(gr==1){
      axis_lbl = sprintf("#%d: No cross-mapping", gr)
    } else{
      q1 = ifelse(gr==2, 0, all_data[[lbl]]$quantile_group_data$quantiles[gr-2])
      q2 = all_data[[lbl]]$quantile_group_data$quantiles[gr-1]
      q_pairs = all_data[[lbl]]$quantile_group_data$n_pairs_quantiles[gr-1]
      axis_lbl = sprintf("#%d: %.4g-%.4g(%d)",gr, q1, q2, q_pairs)
    }
    
    arg_labels[length(arg_labels)+1] = axis_lbl
    arg_colors[length(arg_colors)+1] = color_opts[lbl]
  }
}


names(arg_list) = 'x'
arg_list[['names']] = arg_labels
arg_list[['col']] = arg_colors
arg_list['pchMed'] = '.:'

do.call(colored_vioplot, args = arg_list)
title(main=paste0('Correlation between gene pairs'),
      xlab = '',
      ylab = 'Correlation (r)')
legend("topleft", fill = color_opts,
       legend = data_labels, box.lty=0, bg = rgb(1,1,1,0.5), horiz = T)











wp_mat = matrix(NA, nrow = length(data_labels), ncol = n_groups)
rownames(wp_mat) = data_labels
colnames(wp_mat) = paste0('#', 1:n_groups)
for(lbl in data_labels){
  for(gr in 1:n_groups){
    wp_mat[lbl,gr] = -log10(all_data[[lbl]]$quantile_group_data$quantile_cross_map_wilcox_list[[gr]]$w_r$p.value + 1e-100)
  }
}

barplot(wp_mat, main="Wilcoxon p-values (R_cross > R_non_cross)",
        xlab="Quantile groups", 
        col=color_opts,
        beside=TRUE,
        ylab=expression(-log[10]*p))
legend("topleft", fill = color_opts,
       legend = data_labels, box.lty=0, bg = rgb(1,1,1,0.5), horiz = T)


dev.off()
