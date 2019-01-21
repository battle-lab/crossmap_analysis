# Gencode v19 does not have any sub-category of pseudogenes, but Gencode v26 has.
# But Gencode v26 is based on hg38, rather than hg19.
# We copy the pseudogene category from Gencode 26 to Gencode 19, 
# where the gene id (without version) is same.

library(argparser)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-table", help="table where the pseudogene sub-category is absent", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_analysis/analysis/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt")
args <- add_argument(args, "-annot", help="gene annotation where the pseudogene sub-category is present", default="/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "-genecol", help="column title containing gene id in the table", default="gene")
args <- add_argument(args, "-typecol", help="column title containing gene type in the table", default="gene_type")
args <- add_argument(args, "-o", help="output file", default="results/table_with_pseudogene_subcategory.txt")

argv = parse_args(args)
table_fn = argv$table
annot_fn = argv$annot
gene_col = argv$genecol
gene_type_col = argv$typecol
out_fn = argv$o

### read data
annot_df = read_df(annot_fn, header = T, row.names = F)
table_df = read_df(table_fn, header = T, row.names = F)

### save original column names and order, will require when re-saving output
original_table_cols = colnames(table_df)

### get the gene ids without the version
annot_df$gene_id_wo_version = sapply(annot_df$gene_id, function(g) strsplit(g, split='\\.')[[1]][1])
table_df$gene_id_wo_version = sapply(table_df[,gene_col], function(g) strsplit(g, split='\\.')[[1]][1])

### keep only one row for each gene id (multiple versions are source of multiple rows)
annot_df = annot_df[!duplicated(annot_df$gene_id_wo_version), ]
rownames(annot_df) = annot_df$gene_id_wo_version

is_pseudogenes_in_table = table_df[,gene_type_col] == 'pseudogene'
table_pseudogenes = unique(table_df$gene_id_wo_version[is_pseudogenes_in_table])
pseudo_annot_df = annot_df[annot_df$gene_id_wo_version %in% table_pseudogenes, ]

common_pseudogenes = rownames(pseudo_annot_df)
is_common_pseudogenes_in_table = table_df$gene_id_wo_version %in% common_pseudogenes
table_df[is_common_pseudogenes_in_table, gene_type_col] = pseudo_annot_df[table_df[is_common_pseudogenes_in_table, 'gene_id_wo_version'],'gene_type']

write_df(table_df[,original_table_cols], file= out_fn, row.names = F, col.names = T)

### count pairs from sub-categories of pseudogenes
pseudo_table_df = table_df[is_pseudogenes_in_table, ]
pseudo_table_df[pseudo_table_df$gene_type=='pseudogene', 'gene_type'] = 'not_annotated'
pseudo_counts = sort(table(pseudo_table_df$gene_type), decreasing = T)
pseudo_counts_cross = table(pseudo_table_df$gene_type[pseudo_table_df$cross_mappable_genes != ""])
pseudo_counts_nocross = table(pseudo_table_df$gene_type[pseudo_table_df$cross_mappable_genes == ""])

## plot
plt_fn = gsub(pattern = '.txt$', replacement = '.pseudogene_composition.pdf', x = out_fn)
pdf(plt_fn, width = 6, height = 4, pointsize = 10)
plt_mat = matrix(0, nrow = length(pseudo_counts), ncol=2, dimnames = list(names(pseudo_counts), c('Not cross-mappable', 'Cross-mappable')))
plt_mat[names(pseudo_counts_nocross), 'Not cross-mappable'] = as.numeric(pseudo_counts_nocross)
plt_mat[names(pseudo_counts_cross), 'Cross-mappable'] = as.numeric(pseudo_counts_cross)

color_opts = paste0('gray', 100 - 1:length(pseudo_counts) * 7)
names(color_opts) = names(pseudo_counts)
bp = barplot(plt_mat, 
             beside = T, 
             col = color_opts,
             mar = c(0.2, 0.2, 0.2, 1), 
             ylab = 'Number of eQTLs',
             main = 'Classification of eQTLs corresponding to pseudogenes',
             names.arg = c('Not Cross-mappable', 'Cross-mappable'))
legend('topright', legend = names(color_opts), fill = color_opts, bg = rgb(1,1,1,0.5))


max_cat = 5
bp = barplot(plt_mat[1:max_cat,], 
             beside = T, 
             col = color_opts[1:max_cat],
             mar = c(0.2, 0.2, 0.2, 1), 
             #ylim = c(0,1), 
             ylab = 'Number of eQTLs',
             main = 'Classification of eQTLs corresponding to pseudogenes',
             names.arg = c('Not Cross-mappable', 'Cross-mappable'))
legend('topright', legend = names(color_opts[1:max_cat]), fill = color_opts[1:max_cat], bg = rgb(1,1,1,0.5))

dev.off()

