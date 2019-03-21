#######################################################################
### this script prepares GTEx RSEM data for cross-mappaing analysis ###
#######################################################################

### compute gene-level TPMs from RSEM
tissues=('Muscle_Skeletal' 'Skin_Sun_Exposed_Lower_leg' 'Testis' 'Thyroid' 'Whole_Blood'  )
tissue_files=('Muscle-Skeletal' 'Skin-SunExposed(Lowerleg)' 'Testis' 'Thyroid' 'WholeBlood')
rsem_transcript_tpm_dir="/work-zfs/abattle4/lab_data/GTEx_v7/rna_seq_from_portal/transcript_tpm"
rsem_transcript_count_dir="/work-zfs/abattle4/lab_data/GTEx_v7/rna_seq_from_portal/transcript_count"
annot_fn="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.transcript.txt"
prev_normalized_expr_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression"
prev_cov_dir="/work-zfs/abattle4/lab_data/GTEx_v7_cisEQTL/GTEx_Analysis_v7_eQTL_covariates"

rsem_gene_tpm_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/raw"
rsem_normalized_expr_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/processed_expression"
rsem_cov_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/covariates"


cd "/work-zfs/abattle4/ashis/prog/gtex_exploration"
if [ ! -d $rsem_gene_tpm_dir ]; then mkdir -p $rsem_gene_tpm_dir; fi
for(( i=0; i<${#tissues[@]}; i++)) 
do
  tis="${tissues[$i]}"
  tis_fn="${tissue_files[$i]}"
  tpm_fn="${rsem_transcript_tpm_dir}/${tis_fn}.txt"
  count_fn="${rsem_transcript_count_dir}/${tis_fn}.txt"
  gene_tpm_fn="${rsem_gene_tpm_dir}/${tis}.tpm.txt"
  gene_count_fn="${rsem_gene_tpm_dir}/${tis}.count.txt"
  Rscript src/data_process/summarize_tx2gene_rsem.R -tpm "$tpm_fn" \
                                                    -count "$count_fn" \
                                                    -annot "$annot_fn" \
                                                    -otpm "$gene_tpm_fn" \
                                                    -ocount "$gene_count_fn" &
done

### normalize rsem expression
cd "/work-zfs/abattle4/ashis/prog/crossmap_analysis"
if [ ! -d $rsem_normalized_expr_dir ]; then mkdir $rsem_normalized_expr_dir ; fi
for(( i=0; i<${#tissues[@]}; i++)) 
do
  tis="${tissues[$i]}"
  gene_tpm_fn="${rsem_gene_tpm_dir}/${tis}.tpm.txt"
  prev_processed_expr_fn="${prev_normalized_expr_dir}/${tis}.v7.normalized_expression.txt"
  processed_expr_fn="${rsem_normalized_expr_dir}/${tis}.v7.normalized_expression.txt"
  
  Rscript process_gtex_rsem_data.R -tpm "$gene_tpm_fn" \
                                        -prev "$prev_processed_expr_fn" \
                                        -o "$processed_expr_fn" &
done


### compute peer factors
cd "/work-zfs/abattle4/ashis/prog/crossmap_analysis"
if [ ! -d "$rsem_cov_dir" ]; then mkdir "$rsem_cov_dir" ; fi
for(( i=0; i<${#tissues[@]}; i++)) 
do
  tis="${tissues[$i]}"
  processed_expr_fn="${rsem_normalized_expr_dir}/${tis}.v7.normalized_expression.txt"
  prev_cov_fn="${prev_cov_dir}/${tis}.v7.covariates.txt"
  rsem_cov_fn="${rsem_cov_dir}/${tis}.v7.covariates.txt"
  
  Rscript process_gtex_rsem_cov.R -expr "$processed_expr_fn" \
                                        -prev "$prev_cov_fn" \
                                        -o "$rsem_cov_fn" 2>&1 | tee "$rsem_cov_fn.log"  &
done


### replicate current gtex eqtls using rsem data
cd "/work-zfs/abattle4/ashis/prog/crossmap_analysis"
Rscript process_gtex_rsem_data_for_replication.R

rsem_replication_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/replication_in_rsem"
tissues=('Muscle_Skeletal' 'Skin_Sun_Exposed_Lower_leg' 'Thyroid' 'Whole_Blood')

expr_pfx=""
expr_sfx=".v7.normalized_expression.txt"
cov_pfx=""
cov_sfx=".v7.covariates.txt"
eqtl_pfx=""
eqtl_sfx="_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt"

rsem_cov_names="C1,C2,C3,sex,platform"
for((i=1; i<=60; i++))
do
  rsem_cov_names="$rsem_cov_names,InferredCov${i}"
done

fdr=0.05
rapproach="full_random"
run=3
label="RSEM"

rsem_eqtl_dir="$rsem_replication_dir/rsem_eqtl"
rsem_geno_fn="$rsem_replication_dir/rsem_geno.txt"
rsem_snp_annot_fn="$rsem_replication_dir/rsem_snp_annot.txt"
rsem_gene_annot_fn="$rsem_replication_dir/rsem_gene_annot.txt"

if [ ! -d "$rsem_eqtl_dir" ]; then mkdir "$rsem_eqtl_dir" ; fi
for(( i=0; i<${#tissues[@]}; i++)) 
do
  tis="${tissues[$i]}"
  
  echo "replicating eqtls - $tis"
  
  rsem_expr_fn="${rsem_replication_dir}/rsem_expr/${expr_pfx}${tis}${expr_sfx}"
  rsem_cov_fn="${rsem_replication_dir}/rsem_cov/${cov_pfx}${tis}${cov_sfx}"
  eqtl_fn="${rsem_replication_dir}/rnaseqc_eqtls_crossmap_trans_eqtl_cross_chr/${eqtl_pfx}${tis}${eqtl_sfx}"
  replication_fn="${rsem_eqtl_dir}/${tis}_trans_eqtl_replication_in_rsem_fdr_${fdr}_bestvar.txt"
  
  cd "/work-zfs/abattle4/ashis/prog/misc/gtex_eqtl_v9"
  echo Rscript replicate_gtex_eqtl_in_dgn.R  -eqtl "$eqtl_fn" \
                                              -expr "$rsem_expr_fn" \
                                              -geno "$rsem_geno_fn" \
                                              -cov "$rsem_cov_fn" \
                                              -covname "$rsem_cov_names" \
                                              -snp_annot "$rsem_snp_annot_fn" \
                                              -gene_annot "$rsem_gene_annot_fn" \
                                              -best "TRUE" \
                                              -fdr "$fdr" \
                                              -rapproach "$rapproach" \
                                              -run "$run" \
                                              -lab "$label" \
                                              -o "$replication_fn" 2>&1 | tee "$replication_fn.log"
  
  replication_fn="${rsem_eqtl_dir}/${tis}_trans_eqtl_replication_in_rsem_fdr_${fdr}_allvar.txt"
  echo Rscript replicate_gtex_eqtl_in_dgn.R  -eqtl "$eqtl_fn" \
                                              -expr "$rsem_expr_fn" \
                                              -geno "$rsem_geno_fn" \
                                              -cov "$rsem_cov_fn" \
                                              -covname "$rsem_cov_names" \
                                              -snp_annot "$rsem_snp_annot_fn" \
                                              -gene_annot "$rsem_gene_annot_fn" \
                                              -best "FALSE" \
                                              -fdr "$fdr" \
                                              -rapproach "$rapproach" \
                                              -run "$run" \
                                              -lab "$label" \
                                              -o "$replication_fn" 2>&1 | tee "$replication_fn.log"
                                              
  cd /work-zfs/abattle4/ashis/prog/crossmap_analysis
  
  replication_fn="${rsem_eqtl_dir}/${tis}_trans_eqtl_replication_in_rsem_fdr_${fdr}_bestvar.txt"
  plt_fn=$(echo $replication_fn | sed 's/\.txt$/.txt.crossmap.pdf/g' )
  Rscript plot_rsem_pvalue_by_crossmap.R -eqtl "$replication_fn" -lab "$label" -o "$plt_fn"
  
  replication_fn="${rsem_eqtl_dir}/${tis}_trans_eqtl_replication_in_rsem_fdr_${fdr}_allvar.txt"
  plt_fn=$(echo $replication_fn | sed 's/\.txt$/.txt.crossmap.pdf/g' )
  Rscript plot_rsem_pvalue_by_crossmap.R -eqtl "$replication_fn" -lab "$label" -o "$plt_fn"
done
