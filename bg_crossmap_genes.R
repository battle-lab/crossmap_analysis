### this script computes the background rate of cross-mappable genes in expression data

library(argparser)
library(data.table)
source('io_util.R')

args <- arg_parser('program')
args <- add_argument(args, '-expr', 
                     help='expression files - comma separated',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Muscle_Skeletal.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Skin_Sun_Exposed_Lower_leg.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Testis.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Thyroid.v7.normalized_expression.txt,/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt')
args <- add_argument(args, "-cross", 
                     help="cross-mappability file", 
                     default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength.txt")
args <- add_argument(args, '-label',
                     help='Tissue labels (comma-separated)',
                     default='Muscle - Skeletal,Skin - Sun Exposed,Testis,Thyroid,Whole Blood')
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/background_crossmap.txt')

argv <- parse_args(args)
expr_fn_input = argv$expr
crossmap_fn = argv$cross
labels_input = argv$label
out_fn <- argv$o

### parse inputs
parts = strsplit(expr_fn_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
expr_files = parts[is_valid_parts]

parts = strsplit(labels_input, split = ',')[[1]]
is_valid_parts = sapply(parts, function(s) nchar(s)>0)
tissue_labels = parts[is_valid_parts]

stopifnot(length(expr_files) == length(tissue_labels))

### read cross-map
crossmap_df = fread(input = crossmap_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")

bg_crossmap_rate = sapply(expr_files, function(expr_fn){
  print(expr_fn)
  expr_df = read_df(expr_fn, sep = '\t', row.names = T)
  genes = rownames(expr_df)
  
  ### find cross-mappable gene pairs among genes in tisssue
  g1_genes = as.factor(crossmap_df[,1])
  g1_genes_levels = levels(g1_genes)
  g1_genes_levels_in_tissue = g1_genes_levels %in% genes
  g1_in_tissue = g1_genes_levels_in_tissue[as.integer(g1_genes)]
  
  g2_genes = as.factor(crossmap_df[,2])
  g2_genes_levels = levels(g2_genes)
  g2_genes_levels_in_tissue = g2_genes_levels %in% genes
  g2_in_tissue = g2_genes_levels_in_tissue[as.integer(g2_genes)]
  
  ### compute bg cross-mappable rate (directed)
  crossmap_in_tissue_df = crossmap_df[g1_in_tissue & g2_in_tissue, ]
  n_directed_crossmap_in_tissue = nrow(crossmap_in_tissue_df)
  total_directed_pairs_in_tissue = (length(genes) * (length(genes) - 1))
  directed_crossmpa_rate = n_directed_crossmap_in_tissue / total_directed_pairs_in_tissue
  
  ### compute bg cross-mappable rate (undirected)
  g1 = pmin(crossmap_in_tissue_df[,1], crossmap_in_tissue_df[,2])
  g2 = pmax(crossmap_in_tissue_df[,1], crossmap_in_tissue_df[,2])
  undirected_crossmap_in_tissue_df = data.frame(g1 = g1, g2 = g2)
  undirected_crossmap_in_tissue_df = unique(undirected_crossmap_in_tissue_df)
  n_undirected_crossmap_in_tissue = nrow(undirected_crossmap_in_tissue_df)
  total_undirected_pairs_in_tissue = (length(genes) * (length(genes) - 1))/2
  undirected_crossmpa_rate = n_undirected_crossmap_in_tissue / total_undirected_pairs_in_tissue
  
  return(c(directed_crossmpa_rate=directed_crossmpa_rate, undirected_crossmpa_rate=undirected_crossmpa_rate))
})

bg_crossmap_rate_df = as.data.frame(t(bg_crossmap_rate))
bg_crossmap_rate_df$tissue = tissue_labels
bg_crossmap_rate_df$file = as.character(sapply(expr_files, basename))

write_df(bg_crossmap_rate_df, file = out_fn, row.names = F, col.names = T)
