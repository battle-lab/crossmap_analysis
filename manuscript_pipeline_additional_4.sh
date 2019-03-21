########################################################################
### this script is actually similar to manuscript_pipeline.sh         ##
### to run additional analyses using GTEx RSEM data.                  ##
### Note: some parts of codes have been commented out,                ##
### because we dont need to run those parts for additional analyses.  ##
########################################################################

#############################################################
###
### this script run all cross-mappability analyses.
###
### note-1: data to run this pipeline can be downloaded
### from the link given in the readme file in github.
###
### note-2: if you configure the input/output variables 
### at the begining of this script in the configuration section, 
### everything else should run automatically, unless any 
### memory/disk errors occur.
###
### note-3: if you do not have genotype data, or if you want
### to avoid genotype related analysis (which are generally 
### slow), set run_genotype_analysis=false. otherwise set
### run_genotype_analysis=true
###
### note-4: if you do not have dgn data, or if you want
### to avoid dgn related analysis, please set 
### run_dgn_analysis=false. otherwise set
### run_dgn_analysis=true
###
### note-5: some scripts use slurm jobs. please wait 
### for slurm jobs to complete before moving to next steps.
### also, you may need to edit memory/time to run each job.
### if you do not have slurm, please replace 'sbatch' by 'sh'
### in this script to run jobs locally (and ignore harmless 
### errors).
###
### note-6: try to avoid spaces (other characters) in 
### file paths, that need to be quoted in shell script.
###
#############################################################

##################### configuration section starts here #####################
#### specification
# put the directory of this repository here
crossmap_analysis_prog_dir="/work-zfs/abattle4/ashis/prog/misc/cross_mappability"
# put your directory where all analysis data and results will be stored
manuscript_analysis_dir="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_analysis"

run_genotype_analysis=true  # true to run all analysis, false to avoid genotype related analysis
run_dgn_analysis=false       # true to run dgn analysis, false not to run

rsem_expr_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/processed_expression"
cov_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7_rsem/covariates"
#genotype_vcf_fn="/work-zfs/abattle4/lab_data/GTEx_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz"
processed_genotype_dir="/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process"
tissues=('Muscle_Skeletal' 'Skin_Sun_Exposed_Lower_leg' 'Testis' 'Thyroid' 'Whole_Blood'  )
tissue_labels=('Muscle - Skeletal' 'Skin - Sun Exposed' 'Testis' 'Thyroid' 'Whole Blood')
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)


crossmap_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/hg19_cross_mappability_strength.txt"
gene_annot_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt"
mappability_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/hg19_gene_mappability.txt"
# repeat_mask_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/hg19_RepeatMaskerTrack.txt"
overlap_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/hg19_positional_overlap.txt"
hgnc_gene_family_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/hgnc_gene_family_20180417.txt"

#intermediate_data_dir="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/intermediate_inputs"
intermediate_data_dir="/temp_intermediate/"

dgn_expr_fn="/work-zfs/abattle4/lab_data/dgn/data_used_for_eqtl_study/trans_data.txt"
dgn_geno_fn="/work-zfs/abattle4/ashis/progdata/dgn/genotype/final_completed_genotype.txt"
dgn_snp_annot_fn="/work-zfs/abattle4/ashis/progdata/dgn/genotype/snp_annot.txt"
dgn_cov_fn="/work-zfs/abattle4/lab_data/dgn/covariates/Biological_and_hidden_factors.txt"


pancan_eqtl_fn="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/pancan/trans_eQTLs_all_re.txt"

##################### configuration section ends here #####################

root_in_dir="$manuscript_analysis_dir/data_rsem"
root_out_dir="$manuscript_analysis_dir/analysis_rsem_75mer_36mer_2mismatch"

if [ ! -d $root_in_dir ]; then mkdir -p $root_in_dir; fi
if [ ! -d $root_out_dir ]; then mkdir -p $root_out_dir; fi

processed_expr_dir="$root_in_dir/processed_expression"
corrected_expr_dir="$root_in_dir/corrected_expression"
genotype_processing_dir="$root_in_dir/genotype_process"
filtered_crossmap_dir="$root_in_dir/filtered_crossmap"

crossmap_exploration_dir="$root_out_dir/crossmap_exploration"
gtex_bg_crossmap_rate_dir="$root_out_dir/gtex_bg_crossmap"

random_corr_dir="$root_out_dir/random_corr"
crossmap_in_top_correlated_genes_dir="$root_out_dir/crossmap_in_top_correlated_gene_pairs"

root_eqtl_dir="$root_out_dir/trans_eqtl"
root_matrix_eqtl_dir="$root_eqtl_dir/matrix_eqtls"
inter_chr_trans_eqtl_dir="$root_eqtl_dir/trans_eqtl_cross_chr"
eqtl_gene_mappability_threshold=0.8
inter_chr_map_filtered_trans_eqtl_dir="$root_eqtl_dir/trans_eqtl_cross_chr_mappability_${eqtl_gene_mappability_threshold}"
eqtl_cross_mappability_d=1e6  # 1Mb
inter_chr_map_crossmap_filtered_trans_eqtl_dir="$root_eqtl_dir/trans_eqtl_cross_chr_mappability_${eqtl_gene_mappability_threshold}_crossmap_${eqtl_cross_mappability_d}b"
matrix_eqtl_metadata_dir="$root_eqtl_dir/matrix_eqtls_metadata"
cis_d=1e6
best_cis_snp_dir="$root_eqtl_dir/best_cis_snp_per_gene_${cis_d}"
combined_pos_fn="$genotype_processing_dir/combined_maf_0.05_repeat_masked.012.pos"


crossmap_in_diff_family_fn="$filtered_crossmap_dir/hg19_cross_mappability_strength_diff_family.txt"

trnas_eqtl_dirs=("$inter_chr_trans_eqtl_dir" "$inter_chr_map_filtered_trans_eqtl_dir" )
trans_eqtl_fn_pfx=""
trans_eqtl_fn_sfx=('_cross_chr_trans_eqtl_fdr_0.05.txt' '_cross_chr_map_trans_eqtl_fdr_0.05.txt' )
trans_eqtl_pval_fn_sfx=('_cross_chr_trans_eqtl_all_p_1e-5.txt' '_cross_chr_map_trans_eqtl_all_p_1e-5.txt' )

eqtl_crossmap_dir="$root_out_dir/eqtl_crossmap"
bg_eqtl_crossmap_dir="$root_out_dir/trans_eqtl_bg_rate"
all_genes_bg_eqtl_crossmap_fn="$bg_eqtl_crossmap_dir/all_genes/background_eqtl_crossmap.txt"
coding_genes_bg_eqtl_crossmap_fn="$bg_eqtl_crossmap_dir/protein_coding_genes/background_eqtl_crossmap_protein_coding.txt"
crossmap_filtered_trans_eqtl_dir="$eqtl_crossmap_dir/trans_eqtl_cross_chr_crossmap_1Mb"
random_eqtl_assoc_dir="$root_out_dir/random_eqtl_assoc"

processed_dgn_data_dir="$root_in_dir/dgn_processed"
processed_dgn_expr_fn="$processed_dgn_data_dir/dgn_with_ensembl_genes.txt"
dgn_random_corr_dir="$root_out_dir/random_corr_dgn"
dgn_random_corr_fn="$dgn_random_corr_dir/dgn_random_correlation_trans_data.pdf"
replication_dir="$root_out_dir/eqtl_replication"

processed_pancan_data_dir="$root_in_dir/pancan_processed"
pancan_crossmap_dir="$root_out_dir/eqtl_crossmap_pancan"


#### slurm utility
get_slurm_header()
{
  if [ $# -lt 9 ]; then 
    echo "#not enough arguments in ${FUNCNAME}."; 
    return 1; 
  fi
  
  partition=$1
  workdir=$2
  time=$3
  nodes=$4
  ntasks=$5
  mem=$6
  job_name=$7
  out_fn=$8
  err_fn=$9
  
  cmd_header="#"'!'"/bin/sh"
  cmd_header="$cmd_header\n#SBATCH --partition=${partition}"
  cmd_header="$cmd_header\n#SBATCH --workdir=${workdir}"
  cmd_header="$cmd_header\n#SBATCH --time=${time}"
  cmd_header="$cmd_header\n#SBATCH --nodes=${nodes}"
  cmd_header="$cmd_header\n#SBATCH --ntasks=${ntasks}"
  cmd_header="$cmd_header\n#SBATCH --mem=${mem}"
  cmd_header="$cmd_header\n#SBATCH --job-name=${job_name}"
  cmd_header="$cmd_header\n#SBATCH --output=${out_fn}"
  cmd_header="$cmd_header\n#SBATCH --error=${err_fn}"
  cmd_header="$cmd_header\n"
  
  echo $cmd_header
  return 0
}


### process data in our format from gtex v7 processed data
if [ ! -d $processed_expr_dir ]; then mkdir $processed_expr_dir; fi
for tissue in ${tissues[@]}
do
  echo "copying gtex rsem processed data - $tissue"
  rsem_expr_fn="$rsem_expr_dir/${tissue}.v7.normalized_expression.txt"
  processed_expr_fn="$processed_expr_dir/${tissue}.v7.normalized_expression.txt"
  cp "$rsem_expr_fn" "$processed_expr_fn"
done


### correct expression by removing gtex covariates
if [ ! -d $corrected_expr_dir ]; then mkdir $corrected_expr_dir; fi

cd $crossmap_analysis_prog_dir
for tissue in ${tissues[@]}
do
  echo $tissue
  expr_fn="$processed_expr_dir/$tissue.v7.normalized_expression.txt"
  cov_fn="$cov_dir/$tissue.v7.covariates.txt"
  corrected_expr_fn="$corrected_expr_dir/$tissue.v7.corrected.txt"
  Rscript rm_gtex_cov.R  -expr $expr_fn \
                         -cov $cov_fn \
                         -id 1 \
                         -st 2 \
                         -o $corrected_expr_fn &
done


# ### convert bed file txt file 
# for tissue in ${tissues[@]}
# do
#   echo $tissue
#   
#   bed_fn="$processed_expr_dir/$tissue.v7.normalized_expression.bed"
#   txt_fn="$processed_expr_dir/$tissue.v7.normalized_expression.txt"
#   if [ -f $bed_fn ]; then 
#     # blank title above gene names is required for some i/o
#     cut -f1,3- $bed_fn | awk -F $'\t' '{FS="\t"; OFS="\t"} {if(NR==1) $1=""; print $0}' > $txt_fn
#   fi
#   
#   bed_fn="$corrected_expr_dir/$tissue.v7.corrected.bed"
#   txt_fn="$corrected_expr_dir/$tissue.v7.corrected.txt"
#   if [ -f $bed_fn ]; then
#    cut -f1,3- $bed_fn | awk -F $'\t' '{FS="\t"; OFS="\t"} {if(NR==1) $1=""; print $0}' > $txt_fn
#   fi
#   
# done

# ### explore cross-mappability
# if [ ! -d $crossmap_exploration_dir ]; then mkdir $crossmap_exploration_dir; fi
# dist_plt_fn=$(echo $(basename $crossmap_fn) | sed 's/.txt$/_exploration.pdf/g' )
# Rscript explore_cross_mappability.R -cross "$crossmap_fn" -o "$dist_plt_fn"

# ### compute background cross-mappability rate between two genes using gtex genes
# if [ ! -d $gtex_bg_crossmap_rate_dir ]; then mkdir $gtex_bg_crossmap_rate_dir; fi
# 
# expr_files_str=""
# tissue_labels_str=""
# for((i=0; i<${#tissues[@]}; i++))
# do
#   expr_files_str="$expr_files_str,$processed_expr_dir/${tissues[$i]}.v7.normalized_expression.txt"
#   tissue_labels_str="$tissue_labels_str,${tissue_labels[$i]}"
# done
# 
# gtex_bg_crossmap_rate_fn="$gtex_bg_crossmap_rate_dir/bg_crossmap_in_gtex.txt"
# Rscript bg_crossmap_genes.R -expr "$expr_files_str" \
#                             -cross "$crossmap_fn" \
#                             -label "$tissue_labels_str" \
#                             -o "$gtex_bg_crossmap_rate_fn"
# 
# # plot
# gtex_bg_crossmap_rate_plt_fn="$(echo $gtex_bg_crossmap_rate_fn | sed 's/.txt$/.pdf/g')"
# Rscript plot_bg_crossmap_in_gtex.R -bg "$gtex_bg_crossmap_rate_fn" -o "$gtex_bg_crossmap_rate_plt_fn"

### uncorrected random correlation
if [ ! -d $random_corr_dir ]; then mkdir $random_corr_dir; fi

cd $crossmap_analysis_prog_dir
for tissue in ${tissues[@]}
do
  echo $tissue
  expr_fn="$processed_expr_dir/$tissue.v7.normalized_expression.txt"
  Rscript random_cor_bet_cross_mappable_pairs.R -cross $crossmap_fn \
                                                -expr $expr_fn \
                                                -n 1000 \
                                                -n0 10000 \
                                                -q "" \
                                                -sym FALSE \
                                                -sep '1,2,5,10,20,50,100,200,300,400,500,1e8' \
                                                -o "$random_corr_dir/${tissue}_random_correlation_uncorrected.pdf" &
  # TIME: 13 min
done


### corrected random correlation : gtex cov correction
if [ ! -d $random_corr_dir ]; then mkdir $random_corr_dir; fi

cd $crossmap_analysis_prog_dir
for tissue in ${tissues[@]}
do
  echo $tissue
  expr_fn="$corrected_expr_dir/$tissue.v7.corrected.txt"
  Rscript random_cor_bet_cross_mappable_pairs.R -cross $crossmap_fn \
                                                -expr $expr_fn \
                                                -n 1000 \
                                                -n0 10000 \
                                                -q "" \
                                                -sym FALSE \
                                                -sep '1,2,5,10,20,50,100,200,300,400,500,1e8' \
                                                -o "$random_corr_dir/${tissue}_random_correlation_gtex_corrected.pdf" &
  # TIME: 13 min
done




### corrected random correlation : excluding crossmapping between genes of same family
# create crossmap file considering genes in same family as not-cross-mappable
if [ ! -d $filtered_crossmap_dir ]; then mkdir $filtered_crossmap_dir; fi
cd $crossmap_analysis_prog_dir
Rscript filter_same_family_crossmap.R -cross $crossmap_fn \
                                      -annot $gene_annot_fn \
                                      -family $hgnc_gene_family_fn \
                                      -o $crossmap_in_diff_family_fn
# TIME: ~10 min

# perform random correlation
if [ ! -d $random_corr_dir ]; then mkdir $random_corr_dir; fi

cd $crossmap_analysis_prog_dir
for tissue in ${tissues[@]}
do
  echo $tissue
  expr_fn="$corrected_expr_dir/$tissue.v7.corrected.txt"
  Rscript random_cor_bet_cross_mappable_pairs.R -cross $crossmap_in_diff_family_fn \
                                                -expr $expr_fn \
                                                -n 1000 \
                                                -n0 10000 \
                                                -q "" \
                                                -sym FALSE \
                                                -sep '1,2,5,10,20,50,100,200,300,400,500,1e8' \
                                                -o "$random_corr_dir/${tissue}_random_correlation_gtex_corrected_diff_gene_family.pdf" &
  # TIME: 13 min
done




### cross-mapping fraction in top correlated gene pairs

# specification for slurm computation cluster, you may ignore if you want to run locally
n_threads=1
partition="shared"
nodes=1
ntasks=1
time="0:30:0"
mem="32GB"

script_dir="$crossmap_in_top_correlated_genes_dir/_scripts"
if [ ! -d $crossmap_in_top_correlated_genes_dir ]; then mkdir $crossmap_in_top_correlated_genes_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

for crossmap_threshold in "1e-6" "100"
do
  for tissue in ${tissues[@]}
  do
    echo "===== crossmap in top correlated gene pairs: $tissue - $crossmap_threshold ====="
    
    expr_fn="$root_in_dir/corrected_expression/$tissue.v7.corrected.txt"
    base_pfx="$tissue.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_${crossmap_threshold}"
    
    out_fn="$crossmap_in_top_correlated_genes_dir/$base_pfx.pdf"
    script_fn="$script_dir/$base_pfx.sh"
    
    job_name="t${tissue: 0 : 3}${crossmap_threshold}"
    cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")
    
    cmd_body="module load R"
    cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
    cmd_body="$cmd_body\nRscript crossmap_frac_in_top_correlated_gene_pairs.R \
                -expr $expr_fn \
                -cross $crossmap_fn \
                -mincross $crossmap_threshold \
                -overlap $overlap_fn \
                -o $out_fn"
    cmd_body="$cmd_body\necho DONE"
    cmd_body="$cmd_body\n"
    
    printf "$cmd_header\n" > $script_fn
    printf "$cmd_body" >> $script_fn
    
    #sbatch $script_fn
    sh $script_fn & 
  done
done

###  once the above scripts are completed, create a combined plot
cd $crossmap_analysis_prog_dir 
for crossmap_threshold in "1e-6" "100"
do
  data_files=""
  labels=""
  for tissue in ${tissues[@]}
  do
    data_files="$data_files,$crossmap_in_top_correlated_genes_dir/$tissue.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_${crossmap_threshold}.RData"
    labels="$labels,$tissue"
  done
  combined_plt_fn="$crossmap_in_top_correlated_genes_dir/combined.v7.corrected_crossmap_in_top_correlated_gene_pairs_exc_overlap_threshold_${crossmap_threshold}.pdf"
  Rscript plot_crossmap_in_top_correlated_gene_pairs.R -data $data_files -label $labels -o $combined_plt_fn
done



############ process genotype data ############

# module load bcftools
# module load vcftools
# module load plink
# module load intel/18.0
# module load htslib  # TODO: error loading htslib

if [ ! -d $genotype_processing_dir ]; then mkdir $genotype_processing_dir; fi
if $run_genotype_analysis
then
  for chr in ${chromosomes[@]}
  do
    echo "chr$chr: copying 012 matrix"
    processed_genotype_fn="${processed_genotype_dir}/chr${chr}_maf_0.05_repeat_masked.012.txt"
    chr_genotype_matrix_fn="${genotype_processing_dir}/chr${chr}_maf_0.05_repeat_masked.012.txt"
    cp "$processed_genotype_fn" "$chr_genotype_matrix_fn"
    
    processed_genotype_pos_fn="${processed_genotype_dir}/chr${chr}_maf_0.05_repeat_masked.012.pos"
    chr_genotype_pos_fn="${genotype_processing_dir}/chr${chr}_maf_0.05_repeat_masked.012.pos"
    cp "$processed_genotype_pos_fn" "$chr_genotype_pos_fn"
  done
fi

### combine all snp positions
pfx="chr"
sfx="_maf_0.05_repeat_masked.012.pos"

if $run_genotype_analysis
then
  if [ ! -f $combined_pos_fn ]
  then 
    for chr in ${chromosomes[@]}
    do
      pos_fn="${genotype_processing_dir}/${pfx}${chr}${sfx}"
      cat $pos_fn >> $combined_pos_fn
    done
  fi
fi

############ run trans-eqtls ###########
script_dir="$root_matrix_eqtl_dir/_scripts"
if [ ! -d $root_eqtl_dir ]; then mkdir $root_eqtl_dir; fi
if [ ! -d $root_matrix_eqtl_dir ]; then mkdir $root_matrix_eqtl_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

partition="shared"
nodes=1
ntasks=1
time="16:0:0"
mem="8GB"

max_snp=1000

for tissue in ${tissues[@]}
do
  echo $tissue
  expr_fn="$processed_expr_dir/$tissue.v7.normalized_expression.txt"
  cov_fn="$cov_dir/$tissue.v7.covariates.txt"
  tissue_eqtl_dir="$root_matrix_eqtl_dir/$tissue"
  
  if [ ! -d $tissue_eqtl_dir ]; then mkdir $tissue_eqtl_dir; fi
  
  for chr in ${chromosomes[@]}
  do
    snp_fn="$genotype_processing_dir/chr${chr}_maf_0.05_repeat_masked.012.txt"
    script_fn="$script_dir/${tissue}_run_eqtl_chr${chr}.sh"
    job_name="e${tissue: 0 : 3}${chr}"
    eqtl_prefix="${tissue_eqtl_dir}/meqtl_${tissue}_chr${chr}"
    cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")
  
    cmd_body="module load R"
    cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
    cmd_body="$cmd_body\nRscript trans_eqtl.R -expr $expr_fn -snp $snp_fn -cov $cov_fn -max_snp $max_snp -o $eqtl_prefix"
    cmd_body="$cmd_body\necho DONE"
    cmd_body="$cmd_body\n"
    
    printf "$cmd_header\n" > $script_fn
    printf "$cmd_body" >> $script_fn
    
    if $run_genotype_analysis
    then
      sbatch $script_fn
    fi
    
  done
done


############ compute trans-eqtl fdr : no filter ###########
script_dir="$inter_chr_trans_eqtl_dir/_scripts"
if [ ! -d $inter_chr_trans_eqtl_dir ]; then mkdir $inter_chr_trans_eqtl_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

filter="trans"
fdr_cutoffs='0.05,0.1,0.2'
n_threads=5
partition="shared"
nodes=1
ntasks=$n_threads
time="8:0:0"
mem="20GB"

for tissue in ${tissues[@]}
do
  echo $tissue
  fn_pattern="$root_matrix_eqtl_dir/$tissue/meqtl_${tissue}*.RData"
  tissue_eqtl_prefix="$inter_chr_trans_eqtl_dir/${tissue}_cross_chr_trans_eqtl"

  script_fn="$script_dir/${tissue}_run_trans_fdr.sh"
  job_name="trans${tissue: 0 : 3}"
  cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")

  cmd_body="module load R"
  cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
  cmd_body="$cmd_body\nRscript filter_matrix_eqtl_results.R -files '$fn_pattern' -filter $filter -annot $gene_annot_fn -crossmap $crossmap_fn -thread $n_threads -fdr '$fdr_cutoffs' -o $tissue_eqtl_prefix"
  cmd_body="$cmd_body\necho DONE"
  cmd_body="$cmd_body\n"
  
  printf "$cmd_header\n" > $script_fn
  printf "$cmd_body" >> $script_fn
  
  if $run_genotype_analysis
  then
    sbatch $script_fn
  fi
  
done


############ compute trans-eqtl fdr : gene mappability filter ###########
script_dir="$inter_chr_map_filtered_trans_eqtl_dir/_scripts"
if [ ! -d $inter_chr_map_filtered_trans_eqtl_dir ]; then mkdir $inter_chr_map_filtered_trans_eqtl_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

filter="map_trans"
fdr_cutoffs='0.05,0.1,0.2'
n_threads=5
partition="shared"
nodes=1
ntasks=$n_threads
time="5:0:0"
mem="20GB"

for tissue in ${tissues[@]}
do
  echo $tissue
  fn_pattern="$root_matrix_eqtl_dir/$tissue/meqtl_${tissue}*.RData"
  tissue_eqtl_prefix="$inter_chr_map_filtered_trans_eqtl_dir/${tissue}_cross_chr_${filter}_eqtl"

  script_fn="$script_dir/${tissue}_run_trans_fdr.sh"
  job_name="trans${tissue: 0 : 3}"
  cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")

  cmd_body="module load R"
  cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
  cmd_body="$cmd_body\nRscript filter_matrix_eqtl_results.R -files '$fn_pattern' -filter $filter -annot $gene_annot_fn -map $mappability_fn -crossmap $crossmap_fn -thread $n_threads -min_map $eqtl_gene_mappability_threshold -crossmap_d $eqtl_cross_mappability_d -fdr '$fdr_cutoffs' -o $tissue_eqtl_prefix"
  cmd_body="$cmd_body\necho DONE"
  cmd_body="$cmd_body\n"

  printf "$cmd_header\n" > $script_fn
  printf "$cmd_body" >> $script_fn

  if $run_genotype_analysis
  then
    sbatch $script_fn
  fi
done


### generate matrix eqtl metadata
script_dir="$matrix_eqtl_metadata_dir/_scripts"
if [ ! -d $matrix_eqtl_metadata_dir ]; then mkdir $matrix_eqtl_metadata_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

n_threads=6
partition="shared"
nodes=1
ntasks=$n_threads
time="8:0:0"
mem="40GB"

for tissue in ${tissues[@]}
do
  echo $tissue
  pfx="$root_matrix_eqtl_dir/$tissue/meqtl_${tissue}"
  out_fn="$matrix_eqtl_metadata_dir/meqtl_meta_${tissue}.txt"
  
  script_fn="$script_dir/${tissue}_run_metadata.sh"
  job_name="meta${tissue: 0 : 3}"
  cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")

  cmd_body="module load R"
  cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
  cmd_body="$cmd_body\nRscript generate_eqtl_file_metadata.R -pfx '$pfx' -core $n_threads -o $out_fn"
  cmd_body="$cmd_body\necho DONE"
  cmd_body="$cmd_body\n"
  
  printf "$cmd_header\n" > $script_fn
  printf "$cmd_body" >> $script_fn
  
  if $run_genotype_analysis
  then
    sbatch $script_fn
  fi
done

### generate best cis-snp per gene
script_dir="$best_cis_snp_dir/_scripts"
if [ ! -d $best_cis_snp_dir ]; then mkdir $best_cis_snp_dir; fi
if [ ! -d $script_dir ]; then mkdir $script_dir; fi

n_threads=6
partition="shared"
nodes=1
ntasks=$n_threads
time="8:0:0"
mem="40GB"

for tissue in ${tissues[@]}
do
  echo $tissue
  fn_pattern="$root_matrix_eqtl_dir/$tissue/meqtl_${tissue}*.RData"
  out_fn="$best_cis_snp_dir/best_cis_snp_${tissue}.txt"
  
  script_fn="$script_dir/${tissue}_run_best_cis.sh"
  job_name="cis${tissue: 0 : 3}"
  cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")

  cmd_body="module load R"
  cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
  cmd_body="$cmd_body\nRscript best_cis_snp_per_gene.R -files '$fn_pattern' -annot $gene_annot_fn -d $cis_d -core $n_threads -o $out_fn"
  cmd_body="$cmd_body\necho DONE"
  cmd_body="$cmd_body\n"
  
  printf "$cmd_header\n" > $script_fn
  printf "$cmd_body" >> $script_fn
  
  if $run_genotype_analysis
  then
    sbatch $script_fn
  fi
done

### code to use pre-computed eqtls (copy files in place)
if ! $run_genotype_analysis
then
  if [ ! -d $inter_chr_trans_eqtl_dir ]; then mkdir -p $inter_chr_trans_eqtl_dir; fi
  if [ ! -d $inter_chr_map_filtered_trans_eqtl_dir ]; then mkdir -p $inter_chr_map_filtered_trans_eqtl_dir; fi
  cp "$intermediate_data_dir/trans_eqtl/trans_eqtl_cross_chr"/* $inter_chr_trans_eqtl_dir/
  cp "$intermediate_data_dir/trans_eqtl/trans_eqtl_cross_chr_mappability_0.8"/* $inter_chr_map_filtered_trans_eqtl_dir/
  
  if [ ! -d $bg_eqtl_crossmap_dir ] ; then mkdir -p $bg_eqtl_crossmap_dir; fi
  if [ ! -d "$bg_eqtl_crossmap_dir/all_genes" ] ; then mkdir "$bg_eqtl_crossmap_dir/all_genes"; fi
  if [ ! -d "$bg_eqtl_crossmap_dir/protein_coding_genes" ] ; then mkdir "$bg_eqtl_crossmap_dir/protein_coding_genes"; fi
  cp "$intermediate_data_dir/trans_eqtl_bg_rate/all_genes"/*   "$bg_eqtl_crossmap_dir/all_genes/" 
  cp "$intermediate_data_dir/trans_eqtl_bg_rate/protein_coding_genes"/*   "$bg_eqtl_crossmap_dir/protein_coding_genes/" 
  
  cp "$intermediate_data_dir/genotype_process"/* $genotype_processing_dir/
fi


### combine trans-eqtls
tissues_str=${tissues[0]}
for((i=1; i<${#tissues[@]}; i++))
do
  tissues_str="$tissues_str,${tissues[$i]}"
done

for((ti=0; ti<${#trnas_eqtl_dirs[@]}; ti++))
do
  tdir=${trnas_eqtl_dirs[$ti]}
  sfx=${trans_eqtl_fn_sfx[$ti]}
  pfx=$trans_eqtl_fn_pfx
  
  combined_dir="$tdir/combined_eqtls"
  if [ ! -d $combined_dir ] ; then mkdir $combined_dir; fi
  combined_eqtl_pfx="$combined_dir/${pfx}combined${sfx}"
  combined_eqtl_pfx=$(echo $combined_eqtl_pfx | sed 's/.txt//g')
  
  cd "$crossmap_analysis_prog_dir"
  Rscript combined_trans_eqtl_results.R -dir "$tdir" -tissues "$tissues_str" -pfx "$pfx" -sfx "$sfx" -o "$combined_eqtl_pfx"
done


### annotate cross-mappaing of eqtls
if [ ! -d $eqtl_crossmap_dir ]; then mkdir $eqtl_crossmap_dir; fi

for((ti=0; ti < ${#trnas_eqtl_dirs[@]}; ti++))
do
  tdir=${trnas_eqtl_dirs[$ti]}
  sfx=${trans_eqtl_fn_sfx[$ti]}
  pfx=$trans_eqtl_fn_pfx
  combined_eqtl_fn="$tdir/combined_eqtls/${pfx}combined$(echo $sfx | sed 's/.txt$/.all.unique.txt/g')"
  annotated_combined_eqtl_dir="$eqtl_crossmap_dir/$(basename $tdir)"
  
  if [ ! -f $combined_eqtl_fn ] ; then continue; fi
  if [ ! -d $annotated_combined_eqtl_dir ] ; then mkdir $annotated_combined_eqtl_dir; fi
  
  echo "annotating combined eqtls in $(basename $tdir) ..."
  
  annotated_combined_eqtl_fn="$annotated_combined_eqtl_dir/$(basename $combined_eqtl_fn)"
  annotated_combined_eqtl_fn=$(echo $annotated_combined_eqtl_fn | sed 's/.txt$/.crossmap.txt/g')
  cd $crossmap_analysis_prog_dir
  log_fn=$(echo $annotated_combined_eqtl_fn | sed 's/.txt$/.log/g')
  Rscript annoate_eqtl_crossmap.R -eqtl "$combined_eqtl_fn" \
                                  -cross "$crossmap_fn" \
                                  -d $cis_d \
                                  -permute FALSE \
                                  -genperm FALSE \
                                  -gencode "$gene_annot_fn" \
                                  -o "$annotated_combined_eqtl_fn" \
                                  2>&1 | tee "$log_fn" &
  
  ### annotate individual tissues
  cd $crossmap_analysis_prog_dir
  for((tisi=0; tisi<${#tissues[@]}; tisi++ ))
  do
    tissue=${tissues[$tisi]}
    echo "annotating $tissue in $(basename $tdir) ..."
    eqtl_fn="$tdir/${tissue}${trans_eqtl_pval_fn_sfx[$ti]}"
    
    annotated_eqtl_fn="$annotated_combined_eqtl_dir/${tissue}${trans_eqtl_pval_fn_sfx[$ti]}"
    annotated_eqtl_fn=$(echo $annotated_eqtl_fn | sed 's/.txt$/.crossmap.txt/g')
    log_fn=$(echo $annotated_eqtl_fn | sed 's/.txt$/.log/g')
    Rscript annoate_eqtl_crossmap.R -eqtl "$eqtl_fn" \
                                    -cross "$crossmap_fn" \
                                    -d $cis_d \
                                    -permute FALSE \
                                    -genperm FALSE \
                                    -gencode "$gene_annot_fn" \
                                    -o "$annotated_eqtl_fn" \
                                    2>&1 | tee "$log_fn" &
  done
done

### filter protein-coding genes from trans_eqtls and also from mappabiity-filtered trans_eqtls
cd $crossmap_analysis_prog_dir
for((ti=0; ti < ${#trnas_eqtl_dirs[@]}; ti++))
do
  tdir="$eqtl_crossmap_dir/$(basename ${trnas_eqtl_dirs[$ti]})"
  sfx=$(echo ${trans_eqtl_pval_fn_sfx[$ti]} | sed 's/.txt$/.crossmap.txt/g') 
  coding_eqtl_dir="${tdir}_protein_coding"
  
  if [ ! -d $coding_eqtl_dir ] ; then mkdir -p $coding_eqtl_dir; fi
  
  for tissue in ${tissues[@]}
  do
    echo "filtering eqtl data to contain protein coding genes - $tissue ..."
    in_eqtl_fn="$tdir/${tissue}$sfx"
    out_sfx=$(echo $sfx | sed 's/.crossmap.txt$/_protein_coding.crossmap.txt/g')
    out_eqtl_fn="$coding_eqtl_dir/${tissue}$out_sfx"
    Rscript filter_protein_coding_genes_from_eqtl_results.R -eqtl "$in_eqtl_fn" -gencode "$gene_annot_fn" -o "$out_eqtl_fn" &
  done
done


# ### background cross-mappable eqtl rate
# expr_files_str="$processed_expr_dir/${tissues[0]}.v7.normalized_expression.txt"
# tissue_labels_str="${tissue_labels[0]}"
# for((i=1; i<${#tissues[@]}; i++))
# do
#   expr_files_str="$expr_files_str,$processed_expr_dir/${tissues[$i]}.v7.normalized_expression.txt"
#   tissue_labels_str="$tissue_labels_str,${tissue_labels[$i]}"
# done
# genotype_sfx="_maf_0.05_repeat_masked.012.txt"
# 
# cd $crossmap_analysis_prog_dir
# if  $run_genotype_analysis
# then
#   Rscript bg_crossmap_eqtl_rate.R -expr "$expr_files_str" \
#                                   -geno "$genotype_processing_dir" \
#                                   -genosfx "$genotype_sfx" \
#                                   -cross "$crossmap_fn" \
#                                   -label "$tissue_labels_str" \
#                                   -annot "$gene_annot_fn" \
#                                   -o "$all_genes_bg_eqtl_crossmap_fn" \
#                                   2>&1 | tee "$all_genes_bg_eqtl_crossmap_fn.log" &
#                                   
#   
#   Rscript bg_crossmap_eqtl_rate_protein_coding.R  -expr "$expr_files_str" \
#                                                   -geno "$genotype_processing_dir" \
#                                                   -genosfx "$genotype_sfx" \
#                                                   -cross "$crossmap_fn" \
#                                                   -label "$tissue_labels_str" \
#                                                   -annot "$gene_annot_fn" \
#                                                   -o "$coding_genes_bg_eqtl_crossmap_fn" \
#                                                   2>&1 | tee "$coding_genes_bg_eqtl_crossmap_fn.log" &
#                                 
# fi


### crossmap in top eqtl (without filtering, mappabiity>=0.8, and protein-coding for both)
eqtl_crossmap_dirs=("$eqtl_crossmap_dir/trans_eqtl_cross_chr" \
                    "$eqtl_crossmap_dir/trans_eqtl_cross_chr_mappability_0.8" \
                    "$eqtl_crossmap_dir/trans_eqtl_cross_chr_protein_coding" \
                    "$eqtl_crossmap_dir/trans_eqtl_cross_chr_mappability_0.8_protein_coding")
eqtl_fn_prefixes=("" \
                  "" \
                  "" \
                  "")
eqtl_fn_suffixes=("_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt" \
                  "_cross_chr_map_trans_eqtl_all_p_1e-5.crossmap.txt" \
                  "_cross_chr_trans_eqtl_all_p_1e-5_protein_coding.crossmap.txt" \
                  "_cross_chr_map_trans_eqtl_all_p_1e-5_protein_coding.crossmap.txt")
bg_rate_fns=("" \
            "" \
            "" \
            "") # should use bg rate with filtered genes

for ((run=0; run<${#eqtl_crossmap_dirs[@]}; run++))
do
  crossmap_dir="${eqtl_crossmap_dirs[$run]}"
  eqtl_fn_prefix="${eqtl_fn_prefixes[$run]}"
  eqtl_fn_suffix="${eqtl_fn_suffixes[$run]}"
  bg_rate_fn="${bg_rate_fns[$run]}"
  
  bg_eqtl_crossmap_rate=()
  if [[ ! $bg_rate_fn == "" ]]; then
    # 3rd column contains background eqtl crossmap rate
    bg_eqtl_crossmap_rate=($(cut -f3 $all_genes_bg_eqtl_crossmap_fn | tail -n +2))
  else
    bg_eqtl_crossmap_rate=($(seq -1 -1 -${#tissues[@]}))  # array of negative numbers
  fi
  
  out_dir="$crossmap_dir/crossmap_in_top_eqtl_p_1e-5"
  
  if [ ! -d $out_dir ]; then mkdir $out_dir; fi
  
  cd $crossmap_analysis_prog_dir
  for((ti=0; ti<${#tissues[@]}; ti++ ))
  do
    tissue=${tissues[$ti]}
    bg_rate=${bg_eqtl_crossmap_rate[$ti]}
    echo "$(basename $crossmap_dir): computing crossmap frac in $tissue"
    eqtl_fn="$crossmap_dir/${eqtl_fn_prefix}${tissue}${eqtl_fn_suffix}"
    out_fn="$out_dir/${tissue}_crossmap_in_top_eqtls.pdf"
    Rscript crossmap_frac_in_top_eqtls.R -eqtl "$eqtl_fn" \
                                         -gencode "$gene_annot_fn" \
                                         -rate " $bg_rate" \
                                         -o "$out_fn"
  done
  
  ####### ====== combined plot
  data_files=""
  label_param=""
  for((ti=0; ti<${#tissues[@]}; ti++ ))
  do
    tissue=${tissues[$ti]}
    data_files="$data_files,$out_dir/${tissue}_crossmap_in_top_eqtls.pdf.RData"
    label_param="$label_param,${tissue_labels[$ti]}"
  done
  combined_plt_fn="$out_dir/combined_crossmap_in_top_eqtls.pdf"
  Rscript plot_crossmap_frac_in_top_eqtls.R -data "$data_files" -label "$label_param" -o "$combined_plt_fn"
done


# ########## compute fdr after filtering out cross-mappable eqtls ########## 
# ##### filter crossmappable trans-eqtls
# unfilterd_p_1e_5_sfx="$(echo $trans_eqtl_pval_fn_sfx | sed 's/.txt$/.crossmap.txt/g')"
# crossmap_filtered_pfx=""
# crossmap_filtered_sfx="_cross_chr_crossmap_1Mb_trans_eqtl_all_p_1e-5.crossmap.txt"
# 
# script_dir="$crossmap_filtered_trans_eqtl_dir/_scripts"
# 
# cd "$crossmap_analysis_prog_dir"
# 
# ### slurm settings 
# source 'slurm_util.sh'
# partition="shared"
# nodes=1
# ntasks=1
# time="6:0:0"
# mem="8GB"
# 
# ### create directories
# if [[ ! -d $crossmap_filtered_trans_eqtl_dir ]]; then mkdir $crossmap_filtered_trans_eqtl_dir; fi
# if [[ ! -d $script_dir ]]; then mkdir $script_dir; fi
# 
# ### create and run scripts to filter cross-mappable eqtls
# for tissue in ${tissues[@]}
# do
#   eqtl_fn="$eqtl_crossmap_dir/trans_eqtl_cross_chr/${tissue}$unfilterd_p_1e_5_sfx"
#   expr_fn="$processed_expr_dir/${tissue}.v7.normalized_expression.txt"
#   filtered_eqtl_fn="$crossmap_filtered_trans_eqtl_dir/$crossmap_filtered_pfx${tissue}$crossmap_filtered_sfx"
#   
#   script_fn="$script_dir/filter_crossmap_from_eqtl_results_${tissue}.sh"
#   job_name="fil${tissue: 0 : 3}"
#   cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")
# 
#   cmd_script="module load R\n"
#   cmd_script="${cmd_script}cd \"$crossmap_analysis_prog_dir\"\n"
#   cmd_script="${cmd_script}Rscript filter_crossmap_from_eqtl_results.R -eqtl \"$eqtl_fn\" -gencode \"$gene_annot_fn\" -snp_pos \"$combined_pos_fn\" -expr \"$expr_fn\" -crossmap \"$crossmap_fn\" -d $eqtl_cross_mappability_d -o \"$filtered_eqtl_fn\" \n"
#   cmd_script="${cmd_script}echo DONE\n"
#   
#   printf "$cmd_header" > $script_fn
#   printf "$cmd_script" >> $script_fn
#   
#   #sbatch $script_fn
#   sh $script_fn  &
# done
# 
# # create files with eqtl hits (fdr<=0.05)
# hit_pfx=""
# hit_sfx="_cross_chr_crossmap_1Mb_trans_eqtl_fdr_0.05.crossmap.txt"
# 
# cd $crossmap_analysis_prog_dir
# for tissue in ${tissues[@]}
# do
#   filtered_eqtl_fn="$crossmap_filtered_trans_eqtl_dir/$crossmap_filtered_pfx${tissue}$crossmap_filtered_sfx"
#   hit_eqtl_fn="$crossmap_filtered_trans_eqtl_dir/$hit_pfx${tissue}$hit_sfx"
#   awk 'NR==1 || $5<0.05{print}' $filtered_eqtl_fn > $hit_eqtl_fn
# done
# 
# 
# # combine all eqtl hits at new fdr <=0.05
# tissues_str=$(printf '%s,' ${tissues[@]})
# combined_eqtl_pfx="$crossmap_filtered_trans_eqtl_dir/combined_cross_chr_crossmap_1Mb_trans_eqtl_fdr_0.05"
# 
# Rscript combined_trans_eqtl_results.R -dir "$crossmap_filtered_trans_eqtl_dir" \
#                                       -tissues "$tissues_str" \
#                                       -pfx "$hit_pfx" \
#                                       -sfx "$hit_sfx" \
#                                       -o "$combined_eqtl_pfx"
#                  

############ pseudogene composition  ##############
for((ti=0; ti < ${#trnas_eqtl_dirs[@]}; ti++))
do
  tdir=${trnas_eqtl_dirs[$ti]}
  sfx=$(echo ${trans_eqtl_fn_sfx[$ti]} | sed 's/.txt$/.all.unique.crossmap.txt/g')
  echo "plotting eqtl composition: $(basename $tdir)"
  
  annotated_combined_eqtl_fn="$eqtl_crossmap_dir/$(basename $tdir)/combined${sfx}"
  combined_eqtl_plt_fn=$(echo "$annotated_combined_eqtl_fn" | sed 's/.txt$/.pdf/g')
  cd $crossmap_analysis_prog_dir
  Rscript explore_eqtl_cross_mappability.R -eqtl "$annotated_combined_eqtl_fn" -gencode "$gene_annot_fn" -o "$combined_eqtl_plt_fn"
done


# ############## random eqtl calls  ##############
# symmetric=FALSE
# script_dir="$random_eqtl_assoc_dir/_scripts"
# 
# if [ ! -d $random_eqtl_assoc_dir ]; then mkdir $random_eqtl_assoc_dir; fi
# if [ ! -d $script_dir ]; then mkdir $script_dir; fi
# 
# partition="shared"
# nodes=1
# ntasks=4
# time="10:0:0"
# mem="60GB"
# 
# cd $crossmap_analysis_prog_dir
# for tissue in ${tissues[@]}
# do
#   echo "eqtl association: $tissue ..."
#   expr_fn="$processed_expr_dir/$tissue.v7.normalized_expression.txt"
#   meqtl_dir="$root_matrix_eqtl_dir/$tissue"
#   meta_fn="$matrix_eqtl_metadata_dir/meqtl_meta_$tissue.txt"
#   cis_fn="$best_cis_snp_dir/best_cis_snp_$tissue.txt"
#   
#   script_fn="$script_dir/${tissue}_run_random_eqtl_assoc.sh"
#   job_name="rqtl${tissue: 0 : 3}"
#   cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")
# 
#   cmd_body="module load R"
#   cmd_body="$cmd_body\ncd $crossmap_analysis_prog_dir"
#   cmd_body="$cmd_body\nRscript random_eqtl_association_bet_cross_mappable_pairs.R \\
#     -meqtl $meqtl_dir  \\
#     -meta $meta_fn \\
#     -cis $cis_fn \\
#     -cross $crossmap_fn \\
#     -expr $expr_fn \\
#     -map $mappability_fn \\
#     -annot $gene_annot_fn\\
#     -sym $symmetric \\
#     -n $N_rand \\
#     -q $top_n_pair  \\
#     -core $ntasks  \\
#     -o $random_eqtl_assoc_dir/${tissue}_random_eqtl_${cis_d}.pdf"
#   
#   cmd_body="$cmd_body\necho DONE"
#   cmd_body="$cmd_body\n"
#   
#   printf "$cmd_header\n" > $script_fn
#   printf "$cmd_body" >> $script_fn
#   
#   if $run_genotype_analysis
#   then
#     sbatch $script_fn
#   fi
# done
# 
# # plot -- run after finishing the above scripts
# rdata_files=""
# tissue_labels=""
# for tissue in ${tissues[@]}
# do
#   rdata_fn="$random_eqtl_assoc_dir/${tissue}_random_eqtl_${cis_d}.pdf.RData"
#   if [ -f $rdata_fn ]; then
#     rdata_files="$rdata_files,$rdata_fn"
#     tissue_labels="$tissue_labels,$tissue"
#   fi
# done
# plt_fn="$random_eqtl_assoc_dir/combined_random_eqtl_assoc.pdf"
# if $run_genotype_analysis
# then
#   Rscript plot_random_eqtl_association.R -assoc "$rdata_files" -label "$tissue_labels" -o "$plt_fn"
# fi
# 
# ########## process dgn expression #############
# if [ ! -d $processed_dgn_data_dir ]; then mkdir $processed_dgn_data_dir; fi
# cd "$crossmap_analysis_prog_dir"
# if $run_dgn_analysis
# then
#   Rscript process_dgn_data.R -expr "$dgn_expr_fn" \
#                              -gene_annot "$gene_annot_fn" \
#                              -o "$processed_dgn_expr_fn"
# fi
# ########## co-expression in dgn ########## 
# if [ ! -d $dgn_random_corr_dir ]; then mkdir $dgn_random_corr_dir; fi
# 
# if $run_dgn_analysis
# then
#   Rscript random_cor_bet_cross_mappable_pairs.R -cross "$crossmap_fn" \
#                                                 -expr "$processed_dgn_expr_fn" \
#                                                 -map "$mappability_fn" \
#                                                 -sym FALSE \
#                                                 -n0 10000 \
#                                                 -n 1000 \
#                                                 -q "" \
#                                                 -sep "1,2,5,10,20,50,100,200,300,400,500,1e8" \
#                                                 -o  "$dgn_random_corr_fn" \
#                                                 2>&1 | tee "$dgn_random_corr_fn.log"
# fi
# 
# ########## eQTL replication in dgn ########## 
# gtex_blood_eqtl_fn="$eqtl_crossmap_dir/trans_eqtl_cross_chr/Whole_Blood_cross_chr_trans_eqtl_all_p_1e-5.crossmap.txt"
# dgn_replication_fn="$replication_dir/trans_eqtl_cross_chr/eqtl_replication_dgn.txt"
# dgn_log_fn="$replication_dir/trans_eqtl_cross_chr/eqtl_replication_dgn.log"
# 
# if [ ! -d "$replication_dir" ]; then mkdir "$replication_dir"; fi
# if [ ! -d "$replication_dir/trans_eqtl_cross_chr" ]; then mkdir "$replication_dir/trans_eqtl_cross_chr"; fi
# 
# cd $crossmap_analysis_prog_dir
# if $run_genotype_analysis
# then
#   if $run_dgn_analysis
#   then
#     Rscript replicate_eqtls_in_dgn.R -eqtl "$gtex_blood_eqtl_fn" \
#                                      -expr "$dgn_expr_fn" \
#                                      -geno "$dgn_geno_fn" \
#                                      -cov "$dgn_cov_fn" \
#                                      -snp_annot "$dgn_snp_annot_fn" \
#                                      -gene_annot "$gene_annot_fn" \
#                                      -o "$dgn_replication_fn" \
#                                      2>&1 | tee "$dgn_log_fn"
#   fi
# fi
# 
# ########## replication in gtex ########## 
# # gather eqtl stats
# uniq_eqtl_fn="$eqtl_crossmap_dir/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt"
# tissues_str=$(IFS=, eval 'echo "${tissues[*]}"')
# labels_str=$(IFS=, eval 'echo "${tissue_labels[*]}"')
# meqtl_file_format="$root_matrix_eqtl_dir/#TISSUE#/#FILE#"
# meqtlmeta_file_format="$matrix_eqtl_metadata_dir/meqtl_meta_#TISSUE#.txt"
# eqtl_per_tissue_stats_fn="$replication_dir/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.per_tissue_stats.txt"
# 
# gtex_replication_fn="$replication_dir/trans_eqtl_cross_chr/eqtl_replication_gtex.pdf"
# gtex_log_fn="$replication_dir/trans_eqtl_cross_chr/eqtl_replication_gtex.log"
# if [ ! -d "$replication_dir" ]; then mkdir "$replication_dir"; fi
# if [ ! -d "$replication_dir/trans_eqtl_cross_chr" ]; then mkdir "$replication_dir/trans_eqtl_cross_chr"; fi
# 
# cd $crossmap_analysis_prog_dir
# if $run_genotype_analysis
# then
#   Rscript gather_eqtl_stats_in_gtex.R -eqtl $uniq_eqtl_fn \
#                                       -tis $tissues_str \
#                                       -meqtl $meqtl_file_format \
#                                       -meqtlmeta $meqtlmeta_file_format \
#                                       -o $eqtl_per_tissue_stats_fn
# fi
# 
# Rscript replicate_eqtls_in_gtex.R -eqtl "$eqtl_per_tissue_stats_fn" \
#                                   -tis "$tissues_str" \
#                                   -lab "$labels_str" \
#                                   -o "$gtex_replication_fn" \
#                                   2>&1 | tee "$gtex_replication_fn.log"
# 
# ### pancan (crossmap fraction in top eQTLs)
# if [ ! -d "$processed_pancan_data_dir" ]; then mkdir "$processed_pancan_data_dir"; fi
# if [ ! -d "$pancan_crossmap_dir" ]; then mkdir "$pancan_crossmap_dir"; fi
# 
# cd $crossmap_analysis_prog_dir
# Rscript process_pancan_data.R -annot "$gene_annot_fn" \
#                               -pancan "$pancan_eqtl_fn" \
#                               -o "$processed_pancan_data_dir" \
#                               2>&1 | tee "$processed_pancan_data_dir/process_pancan_data.log"
# 
# combined_eqtl_fn="$processed_pancan_data_dir/Pancan_combined_unique.txt"
# annotated_combined_eqtl_fn=$(echo $combined_eqtl_fn | sed 's/.txt$/.crosssmap.txt/g')
# Rscript annoate_eqtl_crossmap.R -eqtl "$combined_eqtl_fn" \
#                                   -cross "$crossmap_fn" \
#                                   -d $cis_d \
#                                   -permute FALSE \
#                                   -genperm FALSE \
#                                   -gencode "$gene_annot_fn" \
#                                   -o "$annotated_combined_eqtl_fn" \
#                                   2>&1 | tee "$annotated_combined_eqtl_fn.log"
#   
# pancan_crossmap_in_top_eqtl_fn="$pancan_crossmap_dir/Pancan_combined_unique_crossmap_in_top_eqtls.pdf"
# Rscript crossmap_frac_in_top_eqtls.R -eqtl "$annotated_combined_eqtl_fn" \
#                                          -gencode "$gene_annot_fn" \
#                                          -rate " -1" \
#                                          -o "$pancan_crossmap_in_top_eqtl_fn"
#                                          
