get_slurm_header()
{
  if [ $# -lt 9 ]; then 
         echo "#not enough arguments in ${FUNCNAME}."; 
       return -1; 
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

