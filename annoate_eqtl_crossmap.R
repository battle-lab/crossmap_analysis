# this script annotates cross-mappability of eqtl (S-G, snp-gene) pairs
# i.e. finds genes nearby S cross-mappable with G.

library(argparser)
library(data.table)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7_symmetric_crossmap/trans_eqtl/trans_eqtl_cross_chr_mappability_0.8/Whole_Blood_cross_chr_map_trans_eqtl_fdr_0.05.txt")
args <- add_argument(args, "-cross", help="cross-mappability file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/hg19_cross_mappability_strength_symmetric_mean.txt")
args <- add_argument(args, "-d", help="max distance between a gene and the snp", default=1e6)
args <- add_argument(args, "-permute", help="permute genes in the eqtl list", default=FALSE)
args <- add_argument(args, "-genperm", help="permute from gencode, as opposed to permuting among eqtl hits", default=FALSE)
args <- add_argument(args, "-gencode", help="gencode annotation file (txt format)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-o", help="annotated eqtl file", default="results/eqtl_annotated.txt")

argv = parse_args(args)
eqtl_fn = argv$eqtl
cross_mappability_fn = argv$cross
D = argv$d
do_permute = as.logical(argv$permute)
gencode_permute = as.logical(argv$genperm)
gencode_fn = argv$gencode
#biomart_host = argv$mart
out_fn = argv$o

n_in_parts = 10000 # number of eqtls to annotate at a time

### read eqtl file
eqtl_df = read.table(eqtl_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")

### read gencode file
gencode_df = read.table(gencode_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")
rownames(gencode_df) = gencode_df$gene_id
gencode_genes = unique(gencode_df$gene_id)


### permute the genes of eqtls, if instructed
permute_eqtl_df <- function(eqtl_df){
  # permute genes in the eqtl data frame
  # note: if G1 is substitued by G2 in an eQTL, substitue every G1 by G2.
  gene_2_idx = tapply(1:nrow(eqtl_df), eqtl_df$gene, min)
  permuted_gene_2_idx = sample(as.numeric(gene_2_idx), size = length(gene_2_idx), replace = F)
  names(permuted_gene_2_idx) = names(gene_2_idx)
  permuted_idx = permuted_gene_2_idx[eqtl_df$gene]
  eqtl_df$gene = eqtl_df$gene[permuted_idx]
  eqtl_df$gene_chr = eqtl_df$gene_chr[permuted_idx]
  eqtl_df$gene_start = eqtl_df$gene_start[permuted_idx]
  eqtl_df$gene_end = eqtl_df$gene_end[permuted_idx]
  return(eqtl_df)
}

global_permute_eqtl_df <- function(eqtl_df, gencode_df){
  # permute genes in the eqtl data frame
  # unlike permute_eqtl_df(), here a gene can be replaced by any gene from gencode
  # note: if G1 is substitued by G2 in an eQTL, substitue every G1 by G2.
  
  eqtl_genes = unique(eqtl_df$gene)
  permuted_genes = sample(gencode_df$gene_id, size = length(eqtl_genes), replace = F)
  names(permuted_genes) = eqtl_genes
  eqtl_df$gene = permuted_genes[eqtl_df$gene]
  eqtl_df$gene_chr = gencode_df[eqtl_df$gene, 'chr']
  eqtl_df$gene_start = gencode_df[eqtl_df$gene, 'start_pos']
  eqtl_df$gene_end = gencode_df[eqtl_df$gene, 'end_pos']
  return(eqtl_df)
}


if(paste0('', do_permute) == 'TRUE'){
  if(paste0('', gencode_permute) == 'TRUE'){
    eqtl_df <- global_permute_eqtl_df(eqtl_df, gencode_df)
  } else {
    eqtl_df <- permute_eqtl_df(eqtl_df)
  }
}


### read cross-mappable genes
cross_map = fread(input = cross_mappability_fn, sep = '\t', header = F, stringsAsFactors = F, data.table = F, quote = "")
if(ncol(cross_map)==2) {
  # old format: only 2 columns. make symmetric, and fill 3rd column by 1 (equal weight)
  cross_map = data.frame(V1=c(cross_map$V1, cross_map$V2), V2=c(cross_map$V2, cross_map$V1), stringsAsFactors = F)
  cross_map$V3 = 1
}

# S-G is cross-mappable, if reads from another gene g near S map to G
# given G, we need to find g such that crossmap(g->G) > 0
cross_mappability = tapply(cross_map$V1, INDEX = cross_map$V2, FUN = c)



# ### get TSS of every gene
# all_genes = unique(c(gencode_genes, cross_map$V1, cross_map$V2))
# all_genes_ensembl_gene_id = unlist(lapply(all_genes, function(x) strsplit(x, split='[.]')[[1]][1]))
# #all_genes_version = unlist(lapply(all_genes, function(x) strsplit(x, split='[.]')[[1]][2]))
# #ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = biomart_host)
# ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = biomart_host)
# gene_attributes = c("ensembl_gene_id",  "chromosome_name", "transcription_start_site")
# #gene_attributes = c("ensembl_gene_id",  "chromosome_name", "transcript_start")
# tss <- getBM(attributes = gene_attributes,
#              filters = "ensembl_gene_id",
#              values = all_genes_ensembl_gene_id,
#              mart = ensembl)
# # x = colnames(tss)
# # x[3] = "transcription_start_site"
# # colnames(tss) = x
# tss <- aggregate(transcription_start_site ~  ensembl_gene_id + chromosome_name, tss,  min)
# #tss <- aggregate(transcript_start ~  ensembl_gene_id + chromosome_name, tss,  min)
# if(nrow(tss) != length(unique(tss$ensembl_gene_id)))
#   stop('multiple tss found for a single gene!!')
# 
# 
# all_genes_df = data.frame(gene_id=all_genes, ensembl_gene_id=all_genes_ensembl_gene_id, stringsAsFactors = F)
# tss = merge(all_genes_df, tss, by = 'ensembl_gene_id')
# tss$chromosome_name = paste0('chr', tss$chromosome_name)
# rownames(tss) = as.character(tss$gene_id)
# if(nrow(tss) != nrow(all_genes_df))
#   stop('could not find TSS of every gene')


#tss with rownames as gene id, columns: chromosome_name, transcription_start_site
tss_values =  as.integer(apply(gencode_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))

tss = data.frame(chromosome_name = gencode_df[gencode_df$gene_id, 'chr'],
                 transcription_start_site = tss_values, 
                 gene_id = gencode_df$gene_id,
                 row.names = gencode_df$gene_id, 
                 stringsAsFactors = F)



# For each snp-gene pair (S,G)
#   get genes cross-mappable with G
#   filter out cross-mappable genes with distance > D
#   find gene-types of cross-mappable genes
#   find gene-type of G
#   add columns - gene_tss, gene_type cross_mappable_genes cross_tss cross_dist_to_snp cross_type

get_snp_chr_pos <- function(s){
  splitted_id = strsplit(s, split = '_')[[1]]
  if('chr' != paste0(strsplit(splitted_id[1], split = '')[[1]][1:3], collapse = "")){
    splitted_id[1] = paste0('chr', splitted_id[1])
  }
  return(list(chr=splitted_id[1], pos=as.numeric(splitted_id[2])))
}


### function to access tss fast
tss_gene_to_idx <- new.env(hash = T)
tss_genes = rownames(tss)
tmp <- sapply(1:length(tss_genes), function(idx) tss_gene_to_idx[[tss_genes[idx]]] <<- idx )

get_tss_entries <- function(genes, cols=NULL){
  tss_indexes = sapply(genes, function(g) tss_gene_to_idx[[g]])
  if(is.null(cols))
    return(tss[tss_indexes,])
  return(tss[tss_indexes,cols])
}

### function to access gencode fast
gencode_gene_to_idx <- new.env(hash = T)
gencode_genes = rownames(gencode_df)
tmp <- sapply(1:length(gencode_genes), function(idx) gencode_gene_to_idx[[gencode_genes[idx]]] <<- idx )

get_gencode_entries <- function(genes, cols=NULL){
  gencode_indexes = sapply(genes, function(g) gencode_gene_to_idx[[g]])
  if(is.null(cols))
    return(gencode_df[gencode_indexes,])
  return(gencode_df[gencode_indexes,cols])
}

final_annotated_eqtl_df = NULL

n_in_parts = min(n_in_parts, nrow(eqtl_df))
eqtl_part_indexes <- sort(unique(c(0, seq(from = n_in_parts, to = nrow(eqtl_df), by = n_in_parts), nrow(eqtl_df))))
eqtl_df = eqtl_df[with(eqtl_df, order(FDR)), ]

for(pidx in 1:(length(eqtl_part_indexes)-1)){
  start_idx = eqtl_part_indexes[pidx]+1
  end_idx = eqtl_part_indexes[pidx+1]
  print(sprintf("%s: processing %d-%d", Sys.time(), start_idx, end_idx))
  
annotated_eqtl_df <- apply(eqtl_df[start_idx:end_idx,], MARGIN = 1, FUN = function(row){
  S = as.character(row['snps'])
  G = as.character(row['gene'])
  
  # prepare row with extra empty values to return when necessary
  row['gene_tss'] =  get_tss_entries(G, 'transcription_start_site')
  row['gene_type'] = get_gencode_entries(G, 'gene_type')
  row['cross_type'] = NA
  row['cross_mappable_genes'] = NA
  row['cross_tss'] = NA
  row['cross_tss_snp_d'] = NA
  
  row['gene_chr'] = get_gencode_entries(G, 'chr')
  row['gene_start'] = get_gencode_entries(G, 'start_pos')
  row['gene_end'] = get_gencode_entries(G, 'end_pos')
  
  S_chr_pos = get_snp_chr_pos(S)
  row['snps_chr'] = S_chr_pos$chr
  row['snps_pos'] = S_chr_pos$pos
  
  if(!'tissue' %in% names(row))
    row['tissue'] = 'dummy_tissue'

  # get genes cross-mappable with G
  cross_mappable_genes = tryCatch(cross_mappability[[G]], error = function(e) return(NULL))
  cross_mappable_genes = intersect(cross_mappable_genes, rownames(tss))
  if(is.null(cross_mappable_genes) || length(cross_mappable_genes) == 0 ){
    return(row)
  }
  
  # filter out cross-mappable genes with distance > D
  cross_tss =  get_tss_entries(cross_mappable_genes)
  cross_tss$d = abs(cross_tss$transcription_start_site - S_chr_pos$pos)
  #cross_tss$d = abs(cross_tss$transcript_start - S_chr_pos$pos)
  cross_tss = cross_tss[(cross_tss$chromosome_name == S_chr_pos$chr) & cross_tss$d <=D, ]
  if(nrow(cross_tss) == 0){
    return(row) 
  }
  
  # find gene-types of cross-mappable genes
  cross_gene_types = get_gencode_entries(cross_tss$gene_id, 'gene_type')
  
  # add columns - gene_tss, gene_type cross_mappable_genes cross_tss cross_dist_to_snp cross_type, snp_chr, snp_pos
  row['cross_mappable_genes'] = paste(cross_tss$gene_id, collapse = ',')
  row['cross_tss'] = paste(cross_tss$transcription_start_site, collapse = ',')
  row['cross_tss_snp_d'] = paste(cross_tss$d, collapse = ',')
  row['cross_type'] = paste(cross_gene_types, collapse = ',')
  return(row)
})

if(class(annotated_eqtl_df) != 'matrix')
  stop('something went wrong in the above apply() function.')

annotated_eqtl_df = t(annotated_eqtl_df)
part_out_fn = sprintf("%s.part%d.txt", out_fn, pidx)
write.table(annotated_eqtl_df, file = part_out_fn, sep = '\t', row.names = F, col.names = T, quote = F, na = '')

final_annotated_eqtl_df = rbind(final_annotated_eqtl_df, annotated_eqtl_df)
rm(annotated_eqtl_df)
gc(reset = T)
}

write.table(final_annotated_eqtl_df, file = out_fn, sep = '\t', row.names = F, col.names = T, quote = F, na = '')

