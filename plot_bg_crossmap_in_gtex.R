library(argparser)
library(data.table)
source('io_util.R')

args <- arg_parser('program')
args <- add_argument(args, '-bg', 
                     help='Background cross-mapping rate file in gtex',
                     default='/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/explore_cross_mappability/background_crossmap.txt')
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/background_crossmap_gtex.pdf')

argv <- parse_args(args)
bg_crossmap_fn = argv$bg
bg_crossmap_plt_fn <- argv$o

theme_col1 = "deepskyblue"

### plot for bacground cross-mappability rate between genes in gtex
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
