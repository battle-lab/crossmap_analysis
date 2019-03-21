source('io_util.R')

expr_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/processed_expression"
geno_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process"
cov_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/covariates"
annot_fn = "/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt"
eqtl_dir = "/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_analysis/analysis/eqtl_crossmap/trans_eqtl_cross_chr"
fdr_th = 0.05
tissues = c('Muscle_Skeletal', 'Skin_Sun_Exposed_Lower_leg', 'Whole_Blood', 'Testis', 'Thyroid')
chromosomes = 1:22

expr_pfx = ""
expr_sfx = ".v7.normalized_expression.txt"
geno_pfx = "chr"
geno_sfx = "_maf_0.05_repeat_masked.012.txt"
cov_pfx = ""
cov_sfx = ".v7.covariates.txt"
eqtl_pfx = ""
eqtl_sfx = "_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt"

common_covariates = c('C1', 'C2', 'C3', 'sex', 'platform')
peer_factor_pfx = 'InferredCov'

replication_data_dir = '/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/replication_in_rsem'
rsem_genotype_fn = sprintf('%s/rsem_geno.txt', replication_data_dir)
rsem_snp_annot_fn = sprintf('%s/rsem_snp_annot.txt', replication_data_dir)
rsem_gene_annot_fn = sprintf('%s/rsem_gene_annot.txt', replication_data_dir)
rsem_expr_dir = sprintf('%s/rsem_expr', replication_data_dir)
rsem_cov_dir = sprintf('%s/rsem_cov', replication_data_dir)
rnaseqc_eqtl_dir = sprintf('%s/rnaseqc_eqtls_crossmap_trans_eqtl_cross_chr', replication_data_dir)

### combinde genotypes and snp annot
combined_genotype_df = NULL
for(chr in chromosomes){
  print(sprintf('reading genotypes - chr%s', chr))
  geno_fn = sprintf('%s/%s%s%s', geno_dir, geno_pfx, chr, geno_sfx)
  geno_df = read_df(geno_fn)
  
  if(is.null(combined_genotype_df)){
    combined_genotype_df = geno_df
  } else {
    combined_genotype_df = rbind(combined_genotype_df, geno_df[,colnames(combined_genotype_df)])
  }
}

write_df(combined_genotype_df, file = rsem_genotype_fn)

### snp annotation
splits = strsplit(rownames(combined_genotype_df), split="_")
snp_chrs = paste0('chr', sapply(splits, function(x) x[1]))
snp_poss = sapply(splits, function(x) x[2])
snps_loc = data.frame(snp=rownames(combined_genotype_df), chr=snp_chrs, pos = snp_poss, stringsAsFactors = F)
write_df(snps_loc, file = rsem_snp_annot_fn, row.names = F, col.names = T)

### rsem expression file
if(!dir.exists(rsem_expr_dir))
  dir.create(rsem_expr_dir)

for(tissue in tissues){
  print(sprintf('processing expression - %s', tissue))
  expr_fn = sprintf('%s/%s%s%s', expr_dir, expr_pfx, tissue, expr_sfx)
  expr_df = read_df(expr_fn)
  expr_df = t(expr_df)
  rsem_expr_fn = sprintf('%s/%s%s%s', rsem_expr_dir, expr_pfx, tissue, expr_sfx)
  write_df(expr_df, file = rsem_expr_fn)
}

### rsem covariate file
if(!dir.exists(rsem_cov_dir))
  dir.create(rsem_cov_dir)

for(tissue in tissues){
  cov_fn = sprintf('%s/%s%s%s', cov_dir, cov_pfx, tissue, cov_sfx)
  cov_df = as.data.frame(t(read_df(cov_fn)))
  n_uniq = sapply(cov_df, function(x)length(unique(x)))
  # remove constant covariates
  if(any(n_uniq<=1)){ 
    const_covariates = names(n_uniq[n_uniq<=1])
    cov_df = cov_df[,-which(colnames(cov_df) %in% const_covariates), drop=F]
  }
  
  rsem_cov_fn = sprintf('%s/%s%s%s', rsem_cov_dir, cov_pfx, tissue, cov_sfx)
  write_df(cov_df, file = rsem_cov_fn)
}

### rsem gene annotation file
gene_annot_fn = "/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt"
gene_annot_df = read_df(gene_annot_fn, row.names = F)
gene_annot_df$gene_name = gene_annot_df$gene_id
write_df(gene_annot_df, file = rsem_gene_annot_fn, row.names=F, col.names=T)

### eqtl with appropriate column names
if(!dir.exists(rnaseqc_eqtl_dir))
  dir.create(rnaseqc_eqtl_dir)

for(tissue in tissues){
  print(sprintf('processing eqtls - %s', tissue))
  eqtl_fn = sprintf('%s/%s%s%s', eqtl_dir, eqtl_pfx, tissue, eqtl_sfx)
  eqtl_df = read_df(eqtl_fn, row.names = F)
  
  cn = colnames(eqtl_df)
  replaced_names = c(snps = 'rsid', 
                     gene = 'gene.gene_id', 
                     pvalue = 'linreg_struct.pval',
                     gene_tss = 'gene.TSS',
                     gene_chr = 'gene.contig',
                     snps_chr = 'locus.contig',
                     snps_pos = 'locus.position')
  for(i in 1:length(replaced_names))
    cn[cn==names(replaced_names)[i]] = replaced_names[i]
  
  colnames(eqtl_df) = cn
  eqtl_df$gene.gene_name = eqtl_df$gene.gene_id
  
  rnaseqc_eqtl_fn = sprintf('%s/%s%s%s', rnaseqc_eqtl_dir, eqtl_pfx, tissue, eqtl_sfx)
  write_df(eqtl_df, file = rnaseqc_eqtl_fn, row.names = F, col.names = T)
}


