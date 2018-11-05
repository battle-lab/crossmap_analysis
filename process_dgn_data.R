### this script process DGN data
### it change gene-symbols in dgn expression data to uniquely-mapped ensembl ids

library(argparser)
source('io_util.R')

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-expr", help="dgn expr data", default="/work-zfs/abattle4/lab_data/dgn/data_used_for_eqtl_study/trans_data.txt")
args <- add_argument(args, "-gene_annot", help="hg19 gene annotation data (txt)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-o", help="output file", default="results/processed_dgn_expr.txt")

argv = parse_args(args)
dgn_expr_fn = argv$expr
gene_annot_fn = argv$gene_annot
out_fn = argv$o

### read data
dgn_expr_df = read_df(dgn_expr_fn, sep = '\t', header = T, row.names = T, stringsAsFactors = F, quote = "")
gene_annot_df = read_df(gene_annot_fn, sep = '\t', header = T, row.names = F, stringsAsFactors = F, quote = "")

print(sprintf('#samples: %d, #genes: %d', nrow(dgn_expr_df), ncol(dgn_expr_df)))
dgn_genes_df = data.frame(dgn_sym=colnames(dgn_expr_df))

### take symbols uniquely mapped to ensembl ids from gene annotation file
n_ensembl_ids_per_gene = tapply(gene_annot_df$gene_id, INDEX = gene_annot_df$gene_name, FUN = length)
genes_with_unique_symbols = names(n_ensembl_ids_per_gene[n_ensembl_ids_per_gene==1])
unique_sym_gene_annot_df = gene_annot_df[gene_annot_df$gene_name %in% genes_with_unique_symbols,  ,drop=F]

### convert dgn symbols to ensembl gene ids in dgn data
dgn_genes_ensembl_annot_df = merge(unique_sym_gene_annot_df, dgn_genes_df, by.x = 'gene_name', by.y = 'dgn_sym', all = F)
print(sprintf('#genes uniquely mapped to ensembl ids: %d', nrow(dgn_genes_ensembl_annot_df)))
mapped_dgn_expr_df = dgn_expr_df[,dgn_genes_ensembl_annot_df$gene_name, drop=F]
dgn_sym_2_ensembl_id = tapply(dgn_genes_ensembl_annot_df$gene_id, dgn_genes_ensembl_annot_df$gene_name, c)
colnames(mapped_dgn_expr_df) = dgn_sym_2_ensembl_id[colnames(mapped_dgn_expr_df)]

### change dimension to (gene x sample)
mapped_dgn_expr_df = t(mapped_dgn_expr_df)

### save to file
write_df(mapped_dgn_expr_df, file=out_fn, sep='\t')
