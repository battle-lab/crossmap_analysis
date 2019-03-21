suppressMessages(library(scales))
peer_fn_pfx = '/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/peer_selection/n_peer_selection_rsem'

data_fn = sprintf('%s.RData', peer_fn_pfx)
load(data_fn)

color_opts = c('darkolivegreen2', 'orchid', 'lightsalmon1', 'deepskyblue', 'slategray4', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'palegreen', 'bisque4', 'mediumpurple2', 'palevioletred3', 'bisque3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'tan1', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray')
color_opts = color_opts[1:length(n_egenes_per_tissue)]
names(color_opts) = sort(names(n_egenes_per_tissue))


plt_fn = sprintf('%s.plot.pdf', peer_fn_pfx)
pdf(plt_fn)

xmax = max(as.numeric(sapply(n_egenes_per_tissue, function(x) max(as.numeric(names(x)))))) 
xmin = min(as.numeric(sapply(n_egenes_per_tissue, function(x) min(as.numeric(names(x)))))) 
ymax = max(as.numeric(sapply(n_egenes_per_tissue, function(x) max(x)))) 
ymin = min(as.numeric(sapply(n_egenes_per_tissue, function(x) min(x))))

for(ti in 1:length(n_egenes_per_tissue)){
  if(ti == 1){
    plot(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
         y = n_egenes_per_tissue[[ti]], 
         col = color_opts[names(n_egenes_per_tissue)[ti]], 
         pch = 19, 
         lty = 1,
         xlim = c(xmin, xmax), 
         ylim = c(ymin,ymax), 
         xlab = "Number of PEER factors",
         ylab = "Number of eGenes")
  } else{
    points(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
           y = n_egenes_per_tissue[[ti]],
           col=color_opts[names(n_egenes_per_tissue)[ti]], 
           pch=19)
  }
  lines(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
        y = n_egenes_per_tissue[[ti]],
        col=color_opts[names(n_egenes_per_tissue)[ti]], 
        lty=1)
}

legend('topleft', 
       legend = sort(names(n_egenes_per_tissue)),
       pch = 19,
       col = color_opts[sort(names(n_egenes_per_tissue))],
       lty = 1,
       bg = alpha('white', 0.5))

dev.off()
