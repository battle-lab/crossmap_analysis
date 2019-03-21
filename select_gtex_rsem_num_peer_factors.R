suppressMessages(library(MatrixEQTL))
suppressMessages(library(scales))
suppressMessages(source('io_util.R'))

expr_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/processed_expression"
geno_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process"
cov_dir = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/covariates"
annot_fn = "/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt"
chromosomes = c(7,14)
max_genes_per_chr = 1e6 # big number for all genes 
cis_d = 1e6
fdr_th = 0.05
tissues = c('Muscle_Skeletal', 'Skin_Sun_Exposed_Lower_leg', 'Whole_Blood', 'Thyroid', 'Testis')
rand_seed = 101
out_fn = '/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/peer_selection/n_peer_selection_rsem'

expr_pfx = ""
expr_sfx = ".v7.normalized_expression.txt"
geno_pfx = "chr"
geno_sfx = "_maf_0.05_repeat_masked.012.txt"
cov_pfx = ""
cov_sfx = ".v7.covariates.txt"


common_covariates = c('C1', 'C2', 'C3', 'sex', 'platform')
peer_factor_pfx = 'InferredCov'

### read annotation 
annot_df = read_df(annot_fn, row.names = F)
rownames(annot_df) = annot_df[,1]

### function to read covariates, file: cov x sample, output: sample x cov
read_cov <- function(cov_fn, max_factor_levels = 5){
  cov_df = as.data.frame(t(read_df(cov_fn)))
  n_uniq = sapply(cov_df, function(x)length(unique(x)))
  # convert limited option variables to factors
  for(covariate in names(n_uniq[n_uniq<=max_factor_levels])){
    cov_df[,covariate] = as.factor(cov_df[,covariate])
  }
  # remove constant covariates
  if(any(n_uniq<=1)){ 
    const_covariates = names(n_uniq[n_uniq<=1])
    cov_df = cov_df[,-which(colnames(cov_df) %in% const_covariates), drop=F]
  }
  return(cov_df)
}

### function to select n genes among given genes located in a given chr
select_genes <- function(candidate_genes, chr, n, annot_df, seed = NA){
  if(!is.na(seed))
    set.seed(seed)
  chr_str = ifelse(startsWith(as.character(chr), 'chr'), chr, sprintf('chr%s',chr))
  candidate_annot_df = annot_df[candidate_genes,,drop=F]
  chr_candidate_genes = rownames(candidate_annot_df)[candidate_annot_df$chr == chr_str]
  selected_genes = sample(chr_candidate_genes, size = min(n,length(chr_candidate_genes)), replace = F)
  return(selected_genes)
}

n_egenes_per_tissue = lapply(tissues, function(tissue){
  expr_fn = sprintf('%s/%s%s%s', expr_dir, expr_pfx, tissue, expr_sfx)
  expr_df = read_df(expr_fn)
  
  cov_fn = sprintf('%s/%s%s%s', cov_dir, cov_pfx, tissue, cov_sfx)
  cov_df = read_cov(cov_fn)
  #cov_df = cov_df[colnames(expr_df), ]
  
  stopifnot(all(rownames(expr_df) %in% annot_df$gene_id ))
  
  tissue_common_covariates = intersect(common_covariates, colnames(cov_df))
  max_peer_factors = ncol(cov_df) - length(common_covariates)
  n_peer_factors = sort(unique(c(seq(5, max_peer_factors, by = 5), max_peer_factors)))
  n_egenes_per_peer_factors = sapply(n_peer_factors, function(n_peer){
    per_chr_cis_assoc = lapply(chromosomes, function(chr){
      print(sprintf('%s - %d factors - chr%s', tissue, n_peer, chr))
      geno_fn = sprintf('%s/%s%s%s', geno_dir, geno_pfx, chr, geno_sfx)
      snp_df = read_df(geno_fn)
      
      common_samples = intersect(colnames(snp_df), colnames(expr_df))
      common_samples = intersect(rownames(cov_df), common_samples)
      
      selected_genes = select_genes(candidate_genes = rownames(expr_df), chr = chr, n = max_genes_per_chr, annot_df = annot_df, seed = rand_seed + which(tissues == tissue))
      
      nth_peer_covariates = intersect(colnames(cov_df), c(common_covariates, paste0(peer_factor_pfx, 1:n_peer)))
      cov_mat = t(model.matrix( ~ ., data=cov_df[common_samples,nth_peer_covariates]))  # cov x sample
      cov_mat = cov_mat[-which(rownames(cov_mat)=='(Intercept)'),]
      meqtl_cov_slice = SlicedData$new(cov_mat)
      
      meqtl_expr_slice = SlicedData$new(as.matrix(expr_df[selected_genes,common_samples,drop=F]))
      meqtl_snp_slice = SlicedData$new(as.matrix(snp_df[,common_samples]))
      
      gene_loc = annot_df[selected_genes, c('gene_id', 'chr', 'start_pos', 'end_pos')]
      snps_chr_pos = sapply(rownames(snp_df), function(snpid) strsplit(snpid, split='_')[[1]][1:2])
      snps_loc = data.frame(snp=colnames(snps_chr_pos), chr=paste0('chr', snps_chr_pos[1,]), pos = as.numeric(snps_chr_pos[2,]), stringsAsFactors = F)
      
      me = Matrix_eQTL_main(
        snps = meqtl_snp_slice,
        gene = meqtl_expr_slice,
        cvrt = meqtl_cov_slice,
        output_file_name = NULL,
        output_file_name.cis = NULL,
        pvOutputThreshold = 0,
        pvOutputThreshold.cis = 1,
        useModel = modelLINEAR, 
        verbose = FALSE,
        snpspos = snps_loc, 
        genepos = gene_loc,
        cisDist = cis_d,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      meqtl_df = me$cis$eqtls
      # correct pvalue within each gene
      
      best_variant_per_gene = tapply(1:nrow(meqtl_df), INDEX = meqtl_df$gene, function(indexes){
        g_meqtl_df = meqtl_df[indexes, ]
        g_meqtl_df$gene_FDR = p.adjust(g_meqtl_df$pvalue, method = 'BH')
        min_idx = which.min(g_meqtl_df$gene_FDR)
        return(g_meqtl_df[min_idx, ,drop=F])
      })
    
      best_variant_per_gene_df = do.call(rbind, args = best_variant_per_gene)
      return(best_variant_per_gene_df)
    })
    
    all_chr_cis_assoc = do.call(rbind, args = per_chr_cis_assoc)
    all_chr_cis_assoc$overall_FDR = p.adjust(all_chr_cis_assoc$gene_FDR, method = 'BH')
    n_egens = sum(all_chr_cis_assoc$overall_FDR <= fdr_th)
    return(n_egens)
  })
  
  names(n_egenes_per_peer_factors) = n_peer_factors
  return(n_egenes_per_peer_factors)
})
names(n_egenes_per_tissue) = tissues

# save data
data_fn = sprintf('%s.RData', out_fn)
save(n_egenes_per_tissue, file = data_fn)

# plot
plt_fn = sprintf('%s.pdf', out_fn)
pdf(plt_fn)

xmax = max(as.numeric(sapply(n_egenes_per_tissue, function(x) max(as.numeric(names(x)))))) 
xmin = min(as.numeric(sapply(n_egenes_per_tissue, function(x) min(as.numeric(names(x)))))) 
ymax = max(as.numeric(sapply(n_egenes_per_tissue, function(x) max(x)))) 
ymin = min(as.numeric(sapply(n_egenes_per_tissue, function(x) min(x))))

for(ti in 1:length(n_egenes_per_tissue)){
  if(ti == 1){
    plot(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
         y = n_egenes_per_tissue[[ti]], 
         col = ti, 
         pch = ti, 
         lty = ti, 
         xlim = c(xmin, xmax), 
         ylim = c(ymin,ymax), 
         xlab = "Number of PEER factors",
         ylab = "Number of eGenes")
  } else{
    points(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
           y = n_egenes_per_tissue[[ti]],
           col=ti, 
           pch=ti)
  }
  lines(x = as.numeric(names(n_egenes_per_tissue[[ti]])), 
        y = n_egenes_per_tissue[[ti]],
        col=ti, 
        lty=ti)
}

legend('topleft', 
       legend = names(n_egenes_per_tissue),
       pch = 1:length(n_egenes_per_tissue),
       col = 1:length(n_egenes_per_tissue),
       lty = 1:length(n_egenes_per_tissue),
       bg = alpha('white', 0.5))

dev.off()
# 
# for each tissue
#   n_egenes_per_n_peer = NULL
#   for each n_peer
#     n_peer_cis_assoc = NULL
#     for each chr
#       select genes from chr
#       create annotation data for selected genes
#       create expression data for selected genes
#       create genotype data [#genotype_fn: chr${chr}_maf_0.05_repeat_masked.012.txt]
#       call cis-eQTL
#       save to n_peer_cis_assoc
#     
#     correct p-values for number of tests from each gene and take best variant per gene
#     correct p-values for total number of genes
#     compute the number of eGenes
#     save number of eGenes in n_egenes_per_n_peer
# 
# draw plot
# 
# 
