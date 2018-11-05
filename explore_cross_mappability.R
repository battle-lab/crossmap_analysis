library(data.table)
library(argparser)

args <- arg_parser('program')
args <- add_argument(args, '-cross', 
                     help='cross-mappability file',
                     default='/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt')
args <- add_argument(args, '-o',
                     help='output file (pdf)',
                     default='results/explore_cross_mappability_hg19.pdf')

argv <- parse_args(args)
cross_map_fn = argv$cross
plt_fn <- argv$o


### read cross-map data
cross_map = fread(input = cross_map_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
dim(cross_map)
head(cross_map)

summary(cross_map$V3)
summary(cross_map$V4)

### distribution of strength
pdf(plt_fn)
hist(cross_map$V3, breaks = 50, xlab = "Cross-mappability", main = "")
hist(log2(cross_map$V3), breaks = 50, xlab = "log2(Cross-mappability)", col='gray90', main = "")
dev.off()
