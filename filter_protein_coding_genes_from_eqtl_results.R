# filter protein coding genes from eqtl results

library(argparser)
library(data.table)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_mean/trans_eqtl/trans_eqtl_cross_chr/Muscle_Skeletal_cross_chr_trans_eqtl_all_p_1e-5.txt")
args <- add_argument(args, "-gencode", help="gencode annotation file (txt format)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-o", help="annotated eqtl file", default="results/eqtl_protein_coding.txt")

argv = parse_args(args)
eqtl_fn = argv$eqtl
gencode_fn = argv$gencode
out_fn = argv$o

### read eqtl file
eqtl_df = read.table(eqtl_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")

### read gencode file
gencode_df = read.table(gencode_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")

### filter coding genes
coding_genes = gencode_df[gencode_df$gene_type == 'protein_coding', 'gene_id']
coding_eqtl_df = eqtl_df[eqtl_df$gene %in% coding_genes, ]

### write in file
write.table(coding_eqtl_df, file = out_fn, sep = '\t', row.names = F, col.names = T, quote = F)
