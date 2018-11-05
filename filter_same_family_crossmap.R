library(argparser)
library(data.table)

args <- arg_parser('program')
args <- add_argument(args, "-cross", 
                     help="cross-mapping file", 
                     default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt")
args <- add_argument(args, '-annot',
                     help='gene annotation file',
                     default='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, '-family',
                     help='hgnc gene family file',
                     default='/work-zfs/abattle4/lab_data/hgnc_gene_family/hgnc_gene_family_20180417.txt')
args <- add_argument(args, '-o',
                     help='output file.',
                     default='results/diff_family_crossmap.txt')


argv <- parse_args(args)
crossmap_fn = argv$cross
gene_annot_fn = argv$annot
gene_family_fn = argv$family
out_fn <- argv$o

# crossmap_fn = '/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt'
# gene_family_fn='/work-zfs/abattle4/lab_data/hgnc_gene_family/hgnc_gene_family_20180417.txt'
# gene_annot_fn='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt'

cross_map = fread(input = crossmap_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
dim(cross_map)
head(cross_map)

gene_family_df = read.table(gene_family_fn, sep = '\t', header = T, quote = "", comment.char = "", check.names = T)
dim(gene_family_df)
head(gene_family_df)

gene_annot_df = fread(input = gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F, data.table = F, quote = "")
rownames(gene_annot_df) = gene_annot_df$gene_id

### get gene symbols
gene_fac = as.factor(cross_map$V1)
gene_levels = levels(gene_fac)
gene_levels_symbols = gene_annot_df[gene_levels, 'gene_name']
gene1_symbols = gene_levels_symbols[as.integer(gene_fac)]

gene_fac = as.factor(cross_map$V2)
gene_levels = levels(gene_fac)
gene_levels_symbols = gene_annot_df[gene_levels, 'gene_name']
gene2_symbols = gene_levels_symbols[as.integer(gene_fac)]

#cross_map$sym1 = gene1_symbols
#cross_map$sym2 = gene2_symbols

### get family ids
#gene_2_family = aggregate(gene_family_df$Gene.family.ID, by = list(gene_family_df$Approved.Symbol), c)
gene_2_family = tapply(gene_family_df$Gene.family.ID, gene_family_df$Approved.Symbol, c)

sym_fact = as.factor(gene1_symbols)
sym_levels = levels(sym_fact)
sym_levels_family = lapply(sym_levels, function(lev) gene_2_family[[lev]])
sym1_family = sym_levels_family[as.integer(sym_fact)]

sym_fact = as.factor(gene2_symbols)
sym_levels = levels(sym_fact)
sym_levels_family = lapply(sym_levels, function(lev) gene_2_family[[lev]])
sym2_family = sym_levels_family[as.integer(sym_fact)]

### filter same family crossmap
# sym1_family_length = sapply(sym1_family, length)
# sym2_family_length = sapply(sym2_family, length)
# no_fam = sym1_family_length==0 | sym2_family_length==0
t1=Sys.time()
fam_intersect = mapply(intersect, sym1_family, sym2_family)
#fam_intersect = mcmapply(intersect, sym1_family, sym2_family, mc.cores = 4)
t2 = Sys.time()
print(t2 - t1)

has_intersect = sapply(fam_intersect, length) > 0
filtered_cross_map = cross_map[!has_intersect, ]

### write
write.table(filtered_cross_map, file = out_fn, sep = '\t', row.names = F, col.names = F, quote = F)

###
intersected_cross_map = cross_map[has_intersect, ]
write.table(intersected_cross_map, file = paste0(out_fn, '.same_family.txt'), sep = '\t', row.names = F, col.names = F, quote = F)

pdf(paste0(out_fn, '.same_family.pdf'))

hist(intersected_cross_map$V3, breaks = 50, main = 'cross-mappability in same family')
hist(log2(intersected_cross_map$V3), breaks = 50, main = 'cross-mappability in same family', xlab = 'log2(cross-mappability)')
hist(intersected_cross_map$V3[intersected_cross_map$V3>=50], breaks = 50, main = 'cross-mappability (>=50) in same family')
hist(intersected_cross_map$V3[intersected_cross_map$V3>=200], breaks = 50, main = 'cross-mappability (>=200) in same family')

hist(filtered_cross_map$V3, breaks = 50, main = 'cross-mappability in different family')
hist(log2(filtered_cross_map$V3), breaks = 50, main = 'cross-mappability in different family', xlab = 'log2(cross-mappability)')
hist(filtered_cross_map$V3[filtered_cross_map$V3>=50], breaks = 50, main = 'cross-mappability (>=50) in different family')
hist(filtered_cross_map$V3[filtered_cross_map$V3>=200], breaks = 50, main = 'cross-mappability (>=200) in different family')

dev.off()
