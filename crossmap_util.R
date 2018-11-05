library(data.table)

filter_crossmap_by_genes <- function(crossmap_df, incl.genes){
  # crossmap_df : cross-mappability data frame where first two columns are genes
  gene_fac = as.factor(crossmap_df[,1])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene1_included = gene_levels_included[as.integer(gene_fac)]
  
  gene_fac = as.factor(crossmap_df[,2])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene2_included = gene_levels_included[as.integer(gene_fac)]
  
  return(crossmap_df[gene1_included & gene2_included, , drop=F])
}

filter_crossmap_by_crossmappability <- function(crossmap_df, min.crossmap=-Inf, max.crossmap=Inf){
  # 3rd column contains cross-mappability
  if (min.crossmap > -Inf)
    crossmap_df = crossmap_df[crossmap_df[,3]>= min.crossmap, ]
  if (max.crossmap < Inf)
    crossmap_df = crossmap_df[crossmap_df[,3] <= max.crossmap, ]
  
  return(crossmap_df)
}

dataframe_diff <- function(df1, df2){
  return(df1[!duplicated(rbind(df2, df1))[-seq_len(nrow(df2))], ])
}

filter_crossmap_by_overlap <- function(crossmap_df, overlap_df){
  df1 = crossmap_df[,1:2]
  df2 = overlap_df[,1:2]
  colnames(df2) = colnames(df1)
  filtered_df = dataframe_diff(df1, df2)
  filtered_df = merge(crossmap_df, filtered_df)  # add other columns including mappability
  return(filtered_df)
}

directed_to_undirected_crossmap <- function(crossmap_df, FUN=min){
  genes = sort(unique(c(crossmap_df[,1], crossmap_df[,2])))
  genes1 = factor(crossmap_df[,1], levels = genes)
  genes2 = factor(crossmap_df[,2], levels = genes)
  min_idx = pmin(as.integer(genes1), as.integer(genes2))
  max_idx = pmax(as.integer(genes1), as.integer(genes2))
  
  crossmap_dt = data.table(gene1=genes[min_idx], gene2=genes[max_idx], crossmap = crossmap_df[,3])
  symmetric_crossmap_dt = crossmap_dt[ ,.(symmetric_crossmap=FUN(crossmap)), by=.(gene1, gene2)]
  original_colnames = colnames(crossmap_df)
  crossmap_df = as.data.frame(symmetric_crossmap_dt)
  colnames(crossmap_df) = original_colnames
  return(crossmap_df)
}

crossmap_frac_in_genesets <- function(genes, crossmap_df, mappability_df, overlap_df=NULL, crossmap_threshold=1e-6){
  # genes: vector of ensembl ids
  # crossmap_df: cross-mappability, 1st two columns are genes, 3rd column cross-mappability
  # mappability_df: 1st column gene name, 2nd column mappability
  # overlap_df: first two columns are genes, 3rd column is # ofoverlap positions
  
  ### stop for different reasons
  stopifnot(class(genes)=='character')
  stopifnot(class(crossmap_df) == 'data.frame' && is.character(crossmap_df[,1]) && is.character(crossmap_df[,2]) && is.numeric(crossmap_df[,3]))
  stopifnot(class(mappability_df) == 'data.frame' && is.character(mappability_df[,1]) && is.numeric(mappability_df[,2]))
  stopifnot(is.null(overlap_df) || ( class(overlap_df) == 'data.frame' && is.character(overlap_df[,1]) && is.character(overlap_df[,2])))
  
  
  ### remove genes with NA or unmeasures mappability
  rownames(mappability_df) = mappability_df[,1]
  genes = unique(genes)
  unmeasured_genes = setdiff(genes, mappability_df[,1])
  genes = setdiff(genes, unmeasured_genes)
  na_mappable_genes = genes[is.na(mappability_df[genes,2])]
  genes = setdiff(genes, na_mappable_genes)
  
  ### filter crossmap for genes and cross-mappability
  filtered_crossmap_df = filter_crossmap_by_genes(crossmap_df, incl.genes = genes)
  filtered_crossmap_df = filter_crossmap_by_crossmappability(filtered_crossmap_df, min.crossmap = crossmap_threshold)
  filtered_crossmap_df = directed_to_undirected_crossmap(filtered_crossmap_df, FUN = max)
  
  ### filter crossmap and bg gene-pairs for gene-overlap
  n_gene = length(genes)
  total_pairs = (n_gene * (n_gene-1)/2)
  if(!is.null(overlap_df)){
    filtered_overlap_df = filter_crossmap_by_genes(overlap_df, incl.genes = genes)
    filtered_overlap_df = directed_to_undirected_crossmap(filtered_overlap_df[,1:3], FUN = max)
    filtered_crossmap_df = filter_crossmap_by_overlap(filtered_crossmap_df, filtered_overlap_df)
    total_pairs = total_pairs - nrow(filtered_overlap_df)  # exclude overlapped pairs
  }
  
  ### compute crossmap fraction
  crossmap_frac = nrow(filtered_crossmap_df) / total_pairs
  
  return(list(crossmap_frac = crossmap_frac,
              genes = genes,
              unmapped_genes = c(unmeasured_genes, na_mappable_genes) ))
}
