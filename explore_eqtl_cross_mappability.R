# this script explores cross-mappability of eqtls

library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt")
args <- add_argument(args, "-gencode", help="gencode gene annotation file (txt format)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "-o", help="plot file (pdf)", default="results/eqtl_crossmap_explored.pdf")

argv = parse_args(args)
eqtl_fn = argv$eqtl
gencode_fn = argv$gencode
plt_fn = argv$o

LD_BLOCK_SIZE = 10000

gene_types = c("protein_coding", "pseudogene", "antisense", "lincRNA",  "3prime_overlapping_ncrna", 
               "IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", 
               "IG_V_gene", "IG_V_pseudogene", "miRNA", "misc_RNA", "Mt_rRNA", 
               "Mt_tRNA", "polymorphic_pseudogene", "processed_transcript",  "rRNA", "sense_intronic", 
               "sense_overlapping", "snoRNA", "snRNA", "TR_C_gene", "TR_D_gene", 
               "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene")
color_opts = c('palegreen', 'bisque3', 'lightsalmon1', 'deepskyblue', 'tan1', 
               'paleturquoise1', 'lavenderblush3', 'lightsteelblue1', 'orchid', 'cornsilk4', 
               'lightseagreen', 'burlywood3', 'palevioletred1', 'bisque1', 'brown', 
               'darkolivegreen2', 'bisque4', 'mediumpurple2', 'palevioletred3', 'salmon3', 
               'plum1', 'tomato', 'violet', 'wheat4', 'slategray3', 
               'rosybrown1', 'cornsilk4', 'olivedrab3', 'lightslategray')

# additional colors: 'lightcoral', 'pink1', 'cadetblue2', 'azure3', 'lightblue', 'orchid4', 'plum'

stopifnot(length(gene_types) == length(color_opts))
names(color_opts) = gene_types


### read eqtl file
eqtl_df = read.table(eqtl_fn, sep = '\t', header = T, stringsAsFactors = F, quote = "", comment.char = "")
eqtl_df = eqtl_df[with(eqtl_df, order(FDR)),]

# make sure class(cross_tss_snp_d) is character
eqtl_df$cross_tss_snp_d = as.character(eqtl_df$cross_tss_snp_d)
eqtl_df$cross_tss_snp_d[is.na(eqtl_df$cross_tss_snp_d)] = ""

eqtl_df_extended = eqtl_df
eqtl_df_extended$snp_chr = sapply(eqtl_df_extended$snps, function(snpid) strsplit(snpid, split = "_")[[1]][1])
eqtl_df_extended$snp_pos = as.numeric(sapply(eqtl_df_extended$snps, function(snpid) strsplit(snpid, split = "_")[[1]][2]))
eqtl_df_extended$snp_ld_block = floor(eqtl_df_extended$snp_pos / LD_BLOCK_SIZE)

# eqtl with only one snp from a LD block
ld_row_idx_df = aggregate(1:nrow(eqtl_df_extended), by=list(eqtl_df_extended$tissue, eqtl_df_extended$gene, eqtl_df_extended$snp_chr, eqtl_df_extended$snp_ld_block), min)
eqtl_df_ld = eqtl_df_extended[sort(ld_row_idx_df$x),]

# one eQTL from one tissue only
tissue_row_idx_df = aggregate(1:nrow(eqtl_df_ld), by=list(eqtl_df_ld$gene, eqtl_df_ld$snps), min)
eqtl_df_ld_tissue = eqtl_df_ld[sort(tissue_row_idx_df$x),]

### initialize pdf
pdf(plt_fn)

### count how many pairs have cross-mappability within 1Mb of the snp
n_total = nrow(eqtl_df)
n_no_cross = sum(eqtl_df$cross_mappable_genes=="")
n_cross = n_total - n_no_cross

# count #pairs within 100kb, 200kb, ..., 500kb of the snp
is_any_within_d <- function(s, d){
  # s = string of numbers separated by commas
  # d = max distance
  
  prev_opt = options(warn = -1)
  distances = as.numeric(strsplit(s, split = ',')[[1]])
  options(prev_opt)
  return(any(distances[!is.na(distances)] <= d))
}

n_cross_100k = sum(sapply(eqtl_df$cross_tss_snp_d, is_any_within_d, d=100e3))
n_cross_200k = sum(sapply(eqtl_df$cross_tss_snp_d, is_any_within_d, d=200e3))
n_cross_300k = sum(sapply(eqtl_df$cross_tss_snp_d, is_any_within_d, d=300e3))
n_cross_400k = sum(sapply(eqtl_df$cross_tss_snp_d, is_any_within_d, d=400e3))
n_cross_500k = sum(sapply(eqtl_df$cross_tss_snp_d, is_any_within_d, d=500e3))

plot(x = 0:1, y = 0:1, col = 'white', xlab='', ylab='', axes = F)
legend('left', 
       legend = paste0(c('total #pairs: ', 
                         '#pairs with cross-mappable genes within 1Mb of the snp: ', 
                         '#pairs without any cross-mappable genes within 1Mb of the snp: ', 
                         '#pairs with cross-mappable genes within 100kb of the snp: ', 
                         '#pairs with cross-mappable genes within 200kb of the snp: ', 
                         '#pairs with cross-mappable genes within 300kb of the snp: ', 
                         '#pairs with cross-mappable genes within 400kb of the snp: ', 
                         '#pairs with cross-mappable genes within 500kb of the snp: '), 
                       c(n_total, 
                         n_cross, 
                         n_no_cross,
                         n_cross_100k,
                         n_cross_200k,
                         n_cross_300k,
                         n_cross_400k,
                         n_cross_500k) ))

### gene-type among non-cross-mappable and cross-mappable pairs
gene_type_freq_nocross = as.data.frame(table(eqtl_df$gene_type[eqtl_df$cross_mappable_genes==""]))
colnames(gene_type_freq_nocross) = c('gene_type', 'proportion')
gene_type_freq_nocross$proportion = gene_type_freq_nocross$proportion / sum(gene_type_freq_nocross$proportion)

gene_type_freq_cross = as.data.frame(table(eqtl_df$gene_type[eqtl_df$cross_mappable_genes!=""]))
colnames(gene_type_freq_cross) = c('gene_type', 'proportion')
gene_type_freq_cross$proportion = gene_type_freq_cross$proportion / sum(gene_type_freq_cross$proportion)

gene_type_freq = merge(gene_type_freq_nocross, gene_type_freq_cross, by='gene_type', suffixes = c('_no_cross','_cross'), all = T)
gene_type_freq[is.na(gene_type_freq)] = 0


type_vec = levels(gene_type_freq[,1])
col_vec = color_opts[type_vec]
bp = barplot(as.matrix(gene_type_freq[,2:3]), 
             beside = T, 
             col = col_vec, #as.factor( gene_type_freq[,1]), 
             mar = c(0.2, 0.2, 0.2, 1), 
             ylim = c(0,1), 
             ylab = 'Proportion of gene-type',
             main = 'Gene type proportions',
             names.arg = c('eQTL w/o cross-mappability', 'eQTL w/ cross-mappability'))
legend('topright', legend = type_vec, fill = col_vec, bg = rgb(1,1,1,0.5))


### cross-mappability percentage within each gene type
cross_mappability_per_gene_type = tapply(eqtl_df$cross_mappable_genes!="", eqtl_df$gene_type, sum)
cross_mappability_total_per_gene_type = tapply(eqtl_df$cross_mappable_genes!="", eqtl_df$gene_type, length)
cross_mappability_percentage_per_gene_type = cross_mappability_per_gene_type / cross_mappability_total_per_gene_type;
cross_mappability_percentage_per_gene_type = as.data.frame(cross_mappability_percentage_per_gene_type)

cross_mappability_percentage_per_gene_type_plt_arr = cross_mappability_percentage_per_gene_type[,1]
names(cross_mappability_percentage_per_gene_type_plt_arr) = rownames(cross_mappability_percentage_per_gene_type)
prev_mar = par('mar')
par(mar = c(15,6,1,1))
barplot(cross_mappability_percentage_per_gene_type_plt_arr, 
        ylab = 'Fraction of pairs with cross-mappable genes',
        las=2,
        col = color_opts[names(cross_mappability_percentage_per_gene_type_plt_arr)])
par(mar = prev_mar)

# plot showing only the most frequent gene types
n_frequent_gene_types = 4
most_frequent_gene_types = names(cross_mappability_total_per_gene_type[order(-cross_mappability_total_per_gene_type)])[1:min(length(cross_mappability_total_per_gene_type), n_frequent_gene_types)]
cross_mappability_percentage_per_frequent_gene_type = cross_mappability_percentage_per_gene_type[rownames(cross_mappability_percentage_per_gene_type) %in% most_frequent_gene_types, , drop=F]

cross_mappability_percentage_per_frequent_gene_type_plt_arr = cross_mappability_percentage_per_frequent_gene_type[,1]
names(cross_mappability_percentage_per_frequent_gene_type_plt_arr) = rownames(cross_mappability_percentage_per_frequent_gene_type)
prev_mar = par('mar')
par(mar = c(15,6,1,1))
barplot(cross_mappability_percentage_per_frequent_gene_type_plt_arr,
        ylab = 'Proportion of cross-mappable eQTLs',
        las=2, 
        ylim = c(0,1), 
        col = color_opts[names(cross_mappability_percentage_per_frequent_gene_type_plt_arr)])
par(mar = prev_mar)

### only most frequent gene types among non-cross-mappable and cross-mappable eqtls 
gene_type_freq_most_frequent = gene_type_freq[as.character(gene_type_freq$gene_type) %in% most_frequent_gene_types, ]
type_vec_most_frequent = as.character(gene_type_freq_most_frequent[,1])
col_vec_most_frequent = color_opts[type_vec_most_frequent]
bp = barplot(as.matrix(gene_type_freq_most_frequent[,2:3]), 
             beside = T, 
             col = col_vec_most_frequent, #as.factor( gene_type_freq[,1]), 
             mar = c(0.2, 0.2, 0.2, 1), 
             ylim = c(0,1), 
             ylab = 'Proportion of gene type',
             main = 'Gene type proportions',
             names.arg = c('Not Cross-mappable', 'Cross-mappable'))
legend('topright', legend = type_vec_most_frequent, fill = col_vec_most_frequent, bg = rgb(1,1,1,0.5))



### distance distribution of cross-mappable genes, categorized by gene-type
cross_eqtl_df = eqtl_df[eqtl_df$cross_mappable_genes!="",]
cross_tss_snp_d = tapply(cross_eqtl_df$cross_tss_snp_d, cross_eqtl_df$gene_type, paste0, collapse=',' )


for(i in 1:length(cross_tss_snp_d)){
  gene_type = names(cross_tss_snp_d)[i]
  distances = as.numeric(strsplit(cross_tss_snp_d[i], split = ',')[[1]])
  if (length(distances) < 2)
    next
  dens = density(distances)
  print(max(dens$y))
  if (i==1){
    plot(dens, 
         xlim = c(0,1e6), 
         ylim = c(0, 3e-05),
         xlab = "Distance between SNP and cross-mappable gene's TSS",
         col = color_opts[gene_type], 
         lty = i,
         main = 'Distribution of position of cross-mappable genes')
  } else {
    lines(dens, col = color_opts[gene_type], lty = i)
  }
}

legend('topright', legend = names(cross_tss_snp_d), lwd = 1, lty = 1:length(cross_tss_snp_d), col = color_opts[names(cross_tss_snp_d)], bg = rgb(1,1,1,0.5))

# show only pseudogenes -- with counts only
pseudogenes_idx = which(names(cross_tss_snp_d) == 'pseudogene')
for(i in pseudogenes_idx){
  gene_type = names(cross_tss_snp_d)[i]
  distances = as.numeric(strsplit(cross_tss_snp_d[i], split = ',')[[1]])
  if (length(distances) < 2)
    next
  hist(distances, 
       breaks = 50,
       xlab = "Distance between SNP and cross-mappable gene's TSS", 
       main = "Distribution of position of genes cross-mappable to pseudogenes", 
       col = 'gray90')
}


# show position distributions using one SNP from an LD block
cross_eqtl_df_ld = eqtl_df_ld[eqtl_df_ld$cross_mappable_genes!="",]
cross_tss_snp_d_ld = tapply(cross_eqtl_df_ld$cross_tss_snp_d, cross_eqtl_df_ld$gene_type, paste0, collapse=',' )
pseudogenes_idx = which(names(cross_tss_snp_d_ld) == 'pseudogene')
for(i in pseudogenes_idx){
  gene_type = names(cross_tss_snp_d_ld)[i]
  distances = as.numeric(strsplit(cross_tss_snp_d_ld[i], split = ',')[[1]])
  if (length(distances) < 2)
    next
  hist(distances, 
       breaks = 50,
       xlab = "Distance between SNP and cross-mappable gene's TSS", 
       main = "Distribution of position of genes cross-mappable to pseudogenes\n(LD)",
       col = 'gray90')
}

# show position distributions using an eQTL from one tissue only
cross_eqtl_df_ld_tissue = eqtl_df_ld_tissue[eqtl_df_ld_tissue$cross_mappable_genes!="",]
cross_tss_snp_d_ld_tissue = tapply(cross_eqtl_df_ld_tissue$cross_tss_snp_d, cross_eqtl_df_ld_tissue$gene_type, paste0, collapse=',' )
pseudogenes_idx = which(names(cross_tss_snp_d_ld_tissue) == 'pseudogene')
for(i in pseudogenes_idx){
  gene_type = names(cross_tss_snp_d_ld_tissue)[i]
  distances = as.numeric(strsplit(cross_tss_snp_d_ld_tissue[i], split = ',')[[1]])
  if (length(distances) < 2)
    next
  hist(distances, 
       breaks = 50,
       xlab = "Distance between SNP and cross-mappable gene's TSS", 
       main = "Distribution of position of genes cross-mappable to pseudogenes (LD+Tissue)",
       col = 'gray90')
}

### function to select a set of number of top eqtls
select_top_n_eqtls <- function(total_n_eqtls){
  stopifnot(total_n_eqtls > 0)
  if(total_n_eqtls <= 20)
    return(1:total_n_eqtls)
  if(total_n_eqtls <= 100)
    return(sort(unique(c(seq(10, total_n_eqtls, 10), total_n_eqtls))))
  return(sort(unique(c(seq(10, 90, 10), seq(100, total_n_eqtls, 100), total_n_eqtls))))
}

### fraction of genes for each category in top eQTLs
#unique_gene_types = unique(eqtl_df$gene_type)
unique_gene_types = most_frequent_gene_types
first_ns = select_top_n_eqtls(nrow(eqtl_df))
gene_type_fractions = sapply(unique_gene_types, function(gene_type){
  print(gene_type)
  sapply(first_ns, function(first_n){
    frac = sum(eqtl_df$gene_type[1:first_n] == gene_type)/first_n
  })
})

for(i in 1:length(unique_gene_types)){
  gene_type = unique_gene_types[i]
  x = first_ns
  y = gene_type_fractions[, gene_type]
  if (i==1){
    plot(x = x,
         y = y,
         col = color_opts[gene_type],
         lty = i,
         type = 'l',
         xlim = c(0, max(first_ns)),
         ylim = c(0, 1),
         xlab = "Number of top eQTLs",
         ylab = 'Gene type proportion',
         main = 'Gene type proportions in top eQTLs (sorted by FDR)')
  } else {
    lines(x = x,
          y = y, 
          col = color_opts[gene_type],
          lty = i, 
          type = 'l')
  }
  
  points(x = x,
         y = y,
         col = color_opts[gene_type])
}

legend('topright', legend = unique_gene_types, lwd = 1, lty = 1:length(unique_gene_types), col = color_opts[unique_gene_types], bg = rgb(1,1,1,0.5))

### fraction of genes for each category in top eQTLs without cross-mappability
no_cross_eqtl_df = eqtl_df[eqtl_df$cross_mappable_genes=="",]
first_ns = select_top_n_eqtls(nrow(no_cross_eqtl_df))
gene_type_fractions = sapply(most_frequent_gene_types, function(gene_type){
  print(gene_type)
  sapply(first_ns, function(first_n){
    frac = sum(no_cross_eqtl_df$gene_type[1:first_n] == gene_type)/first_n
  })
})

for(i in 1:length(most_frequent_gene_types)){
  gene_type = most_frequent_gene_types[i]
  x = first_ns
  y = gene_type_fractions[, gene_type]
  if (i==1){
    plot(x = x,
         y = y,
         col = color_opts[gene_type],
         lty = i,
         type = 'l',
         xlim = c(0, max(first_ns)),
         ylim = c(0, 1),
         xlab = "Number of top eQTLs",
         ylab = 'Gene type proportion',
         main = 'Gene type proportions in top eQTLs w/o cross-mappability')
  } else {
    lines(x = x,
          y = y, 
          col = color_opts[gene_type],
          lty = i, 
          type = 'l')
  }
  
  points(x = x,
         y = y,
         col = color_opts[gene_type])
}

legend('topright', legend = most_frequent_gene_types, lwd = 1, lty = 1:length(most_frequent_gene_types), col = color_opts[most_frequent_gene_types], bg = rgb(1,1,1,0.5))

### fraction of genes for each category in top eQTLs with cross-mappability
cross_eqtl_df = eqtl_df[eqtl_df$cross_mappable_genes!="",]
first_ns = select_top_n_eqtls(nrow(cross_eqtl_df))
gene_type_fractions = sapply(most_frequent_gene_types, function(gene_type){
  print(gene_type)
  sapply(first_ns, function(first_n){
    frac = sum(cross_eqtl_df$gene_type[1:first_n] == gene_type)/first_n
  })
})

for(i in 1:length(most_frequent_gene_types)){
  gene_type = most_frequent_gene_types[i]
  x = first_ns
  y = gene_type_fractions[, gene_type]
  if (i==1){
    plot(x = x,
         y = y,
         #col = i,
         col = color_opts[gene_type],
         lty = i,
         type = 'l',
         xlim = c(0, max(first_ns)),
         ylim = c(0, 1),
         xlab = "Number of top eQTLs",
         ylab = 'Gene type proportion',
         main = 'Gene type proportions in top eQTLs w/ cross-mappability')
  } else {
    lines(x = x,
          y = y, 
          col = color_opts[gene_type],
          lty = i, 
          type = 'l')
  }
  
  points(x = x,
         y = y,
         col = color_opts[gene_type])
}

legend('topright', legend = most_frequent_gene_types, lwd = 1, lty = 1:length(most_frequent_gene_types), col = color_opts[most_frequent_gene_types], bg = rgb(1,1,1,0.5))


### fraction of cross-mappabile pairs in top eQTLs
first_ns = select_top_n_eqtls(nrow(eqtl_df))
crossmap_fractions = sapply(first_ns, function(first_n){
  frac = sum(eqtl_df$cross_mappable_genes[1:first_n] != "" )/first_n
  return(frac)
})

x = first_ns
y = crossmap_fractions
plot(x = x,
     y = y,
     col = 'coral3',
     lty = 1,
     type = 'l',
     xlim = c(0, max(first_ns)),
     ylim = c(0, 1),
     xlab = "Number of top eQTLs",
     ylab = 'Proportion of cross-mappable eQTL pairs',
     main = 'Cross-mappable eQTL proportions in top eQTLs w/ cross-mappability')

points(x = x,
       y = y,
       col = 'coral3')

dev.off()
