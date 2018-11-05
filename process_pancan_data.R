########################## PancanQTL ######################################
# - PancanQTL (PancanQTL: systematic identification of cis-eQTLs and trans-eQTLs in 33 cancer types, NAR 2018)
# - 231 210 trans eQTL-gene pairs in 33 cancer types at a per-tissue FDR < 0.05, which corresponded
# - to a median P-value < 1.54 × 10−9
# - rna-seq
# - hg19 used for genomic location
# - trans-eQTLs defined as >1Mb from TSS.
# - Downloaded from here: http://bioinfo.life.hust.edu.cn/PancanQTL/download

source('io_util.R')

library(argparser)
library(data.table)

args <- arg_parser("program");
args <- add_argument(args, '-annot',
                     help='hg19 annotation file',
                     default='/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt')
args <- add_argument(args, '-pancan',
                     help='pancan eqtl file',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/published_studies/PancanQTL/trans_eQTLs_all_re.txt')
args <- add_argument(args, '-o',
                     help='output eqtl directory',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/published_studies/processed')

argv = parse_args(args)

# hg19_gene_annot_fn = '/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt'
# eqtl_fn = '/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/published_studies/PancanQTL/trans_eQTLs_all_re.txt'
# processed_eqtl_out_dir = '/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/published_studies/processed'

hg19_gene_annot_fn = argv$annot
eqtl_fn = argv$pancan
processed_eqtl_out_dir = argv$o


### read data
eqtl_df = read_df(eqtl_fn, sep = '\t', header = T, row.names = F) 
colnames(eqtl_df) = c('Tumor', 'Rs_ID', 'SNP_CHR', 'SNP_POS', 'ALLELES', 'GENE', 'GENE_LOC', 'BETA', 'T', 'P')
eqtl_df$GENE_CHR = sapply(eqtl_df$GENE_LOC, function(s) strsplit(s, split = ':')[[1]][1])
dim(eqtl_df) 
# [1] 231210     11

tumors = sort(unique(eqtl_df$Tumor))

for(tumor in tumors){
  print(tumor)
  processed_eqtl_out_fn = paste0(processed_eqtl_out_dir, '/Pancan.', tumor, ".txt")
  
  ### get eqtls in tumor
  trans_eqtl_df = eqtl_df[eqtl_df$Tumor == tumor, ]
  dim(trans_eqtl_df)
  # [1] 984  11
  
  ### get trans-eqtl in our definition on different chr
  trans_eqtl_df = trans_eqtl_df[trans_eqtl_df$SNP_CHR != trans_eqtl_df$GENE_CHR,]
  dim(trans_eqtl_df)
  # [1] 937   7
  
  ### remove duplicated pairs due to different phenotypes in GWAS, and multiple transcripts per gene
  duplicated_pair = duplicated(trans_eqtl_df[,c('Rs_ID', 'GENE')])
  trans_eqtl_df = trans_eqtl_df[!duplicated_pair, ]
  dim(trans_eqtl_df)
  # [1] 937   7
  
  ### get ensembl gene id from gene symbol
  # remove non-uniquely mapped gene symbols
  hg19_gene_annot_df = read_df(hg19_gene_annot_fn, sep = '\t', header = T, row.names = F)
  n_sym_in_annot = table(hg19_gene_annot_df$gene_name)
  non_unique_sym = names(n_sym_in_annot[n_sym_in_annot>1])
  hg19_gene_annot_df = hg19_gene_annot_df[! hg19_gene_annot_df$gene_name %in% non_unique_sym, ]
  annotated_trans_eqtl_df = merge(trans_eqtl_df, hg19_gene_annot_df, by.x = 'GENE', by.y = 'gene_name')
  dim(annotated_trans_eqtl_df)
  # [1] 623  15 ; successfully mapped to an ensembl gene id
  
  ### exclude entries where the chr of ensembl id and the chr of transcript mismatches
  annotated_trans_eqtl_df = annotated_trans_eqtl_df[annotated_trans_eqtl_df$GENE_CHR == annotated_trans_eqtl_df$chr,]
  dim(annotated_trans_eqtl_df)
  # [1] 623  15
  
  ### create data in our format
  #snps    gene    statistic       pvalue  FDR     beta
  formatted_trans_eqtl_df = annotated_trans_eqtl_df
  formatted_trans_eqtl_df$snps = paste(annotated_trans_eqtl_df$SNP_CHR, annotated_trans_eqtl_df$SNP_POS, sep='_')
  formatted_trans_eqtl_df$gene = formatted_trans_eqtl_df$gene_id
  formatted_trans_eqtl_df$statistic = formatted_trans_eqtl_df$T
  formatted_trans_eqtl_df$pvalue = formatted_trans_eqtl_df$P
  formatted_trans_eqtl_df$FDR = formatted_trans_eqtl_df$P # keeping same as p-value, as #tests is not known
  formatted_trans_eqtl_df$beta = formatted_trans_eqtl_df$BETA
  formatted_trans_eqtl_df = formatted_trans_eqtl_df[,c('snps','gene','statistic','pvalue','FDR','beta',
                                                       'Rs_ID', 'ALLELES', 'GENE_LOC', 'GENE_CHR', 'Tumor')]
  
  
  ### save
  write_df(formatted_trans_eqtl_df, file = processed_eqtl_out_fn, row.names = F, col.names = T)
}


### combine all eqtls -- take unique pairs, pvalue = lowest pvalue across all tissues
combined_unique_pancan_trans_eqtl_fn = paste0(processed_eqtl_out_dir, '/Pancan_combined_unique.txt')
trans_eqtl_df_list = lapply(tumors, function(tumor){
  print(tumor)
  processed_eqtl_fn = paste0(processed_eqtl_out_dir, '/Pancan.', tumor, ".txt")
  processed_eqtl_df = read_df(fn = processed_eqtl_fn, header = T, row.names = F)
  return(processed_eqtl_df)
})
combined_trans_eqtl_df = do.call(rbind, trans_eqtl_df_list)
combined_trans_eqtl_df = combined_trans_eqtl_df[with(combined_trans_eqtl_df, order(pvalue)), ]
first_occur = !duplicated(combined_trans_eqtl_df[,c('snps', 'gene')])
unique_combined_trans_eqtl_df = combined_trans_eqtl_df[first_occur, ]
write_df(unique_combined_trans_eqtl_df, file = combined_unique_pancan_trans_eqtl_fn, row.names = F, col.names=T)
