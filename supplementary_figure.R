source('io_util.R')

theme_col1 = "deepskyblue"

### plot cross-mappability distribution
cross_map_fn = '/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt'
cross_map = read_df(fn = cross_map_fn, header = F, row.names = F)
cross_map_plt_fn = "/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/explore_cross_mappability/explore_cross_mappability_hg19_v2.pdf"
pdf(cross_map_plt_fn, pointsize = 10, width = 4, height = 4)
par(mai = c(1.5, 1, 0.1, 0.1)) 
hist(log2(cross_map$V3), breaks = 50, 
     xlab = expression(paste("log"["2"]~"(Cross-mappability)")), 
     ylab = "Number of gene pairs",
     main = "",
     col=theme_col1)
dev.off()


### plot for bacground cross-mappability rate between genes in gtex
bg_crossmap_fn = '/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/explore_cross_mappability/background_crossmap.txt'
bg_crossmap_plt_fn = '/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/explore_cross_mappability/background_crossmap.pdf'

bg_df = read_df(bg_crossmap_fn, header = T, row.names=F)
pdf(bg_crossmap_plt_fn, pointsize = 10, width = 4, height = 4)
fractions = bg_df$directed_crossmpa_rate * 100
names(fractions) = bg_df$tissue
par(mai = c(1.5, 1, 0.1, 0.1)) 
barplot(fractions,
        main="", 
        xlab="",
        ylab = "Background percentage of\ncross-mappable gene pairs (%)",
        col = theme_col1,
        las = 2)

dev.off()

