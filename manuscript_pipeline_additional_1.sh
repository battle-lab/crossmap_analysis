#####################################################################
### this script run additional analysis on cross-mappability data ###
### to estimate contribution of pseduogene sub-categories,        ###
### to estimate contribution of segdups in cross-mappable eqtls,  ###
### etc.                                                          ###
#####################################################################

### settings
annotated_combined_eqtl_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_analysis/analysis/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt"
pseudo_gene_annot_fn="/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt"
pseudogene_annotated_combined_eqtl_fn=$(echo "$annotated_combined_eqtl_fn" | sed 's/.txt$/.pseudo_gencode26.txt/g')
pseudogene_annotated_combined_eqtl_plt_fn=$(echo "$pseudogene_annotated_combined_eqtl_fn" | sed 's/.txt$/.pdf/g')
segdup_fn="/work-zfs/abattle4/lab_data/hg19_tracks/segdup/GRCh37GenomicSuperDup.tab"
gene_annot_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt"
crossmap_analysis_prog_dir="/work-zfs/abattle4/ashis/prog/crossmap_analysis"

### explore pseduogene sub-categories
cd $crossmap_analysis_prog_dir
Rscript annotate_pseudogene_category.R -table "$annotated_combined_eqtl_fn" -annot "$pseudo_gene_annot_fn" -genecol "gene" -typecol "gene_type" -o "$pseudogene_annotated_combined_eqtl_fn"
Rscript explore_eqtl_cross_mappability.R -eqtl "$pseudogene_annotated_combined_eqtl_fn" -o "$pseudogene_annotated_combined_eqtl_plt_fn"


### get eqtls in segdups
annotated_combined_eqtl_in_segdup_fn="$(echo $annotated_combined_eqtl_fn | sed 's/.txt$/.segdup.txt/g')"
annotated_combined_eqtl_in_segdup_plt_fn="$(echo $annotated_combined_eqtl_in_segdup_fn | sed 's/.txt$/.pdf/g')"
cd $crossmap_analysis_prog_dir
Rscript get_eqtls_in_segdup.R -eqtl "$annotated_combined_eqtl_fn" \
                              -segdup "$segdup_fn" \
                              -o "$annotated_combined_eqtl_in_segdup_fn" 
                                
Rscript explore_eqtl_cross_mappability.R -eqtl "$annotated_combined_eqtl_in_segdup_fn" \
                                         -gencode "$gene_annot_fn" \
                                         -o "$annotated_combined_eqtl_in_segdup_plt_fn"


