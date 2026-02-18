#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=run_wf
#SBATCH --output=slurm_out/workflowr_merscope_human_breast_cancer.out
#SBATCH --error=slurm_out/workflowr_merscope_human_breast_cancer.err
#SBATCH --time=06:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=12

################################################################################
echo "------------------------------------------------------------------------"
echo "Job Started on $(date)"
echo "System: $(uname -a)"
echo "CPU: $(lscpu | grep 'Model name')"
echo "Total Memory: $(free -h | awk '/Mem:/ {print $2}')"
echo "Node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "------------------------------------------------------------------------"

################################################################################

msg "======== begin ========"           #
msg 'Working directory ' `pwd`          # current job working directory
msg 'Job run on nodes ' $SLURM_NODELIST # current job assigned nodes
msg 'Job ntasks assign ' $SLURM_TASKS_PER_NODE #
msg 'NCPU per task' $SLURM_CPUS_PER_TASK # Number of CPUs per task
msg 'Total allocated cores ' $SLURM_CPUS_ON_NODE     # calculated total allocated NCPU
msg 'Job ID ' ${SLURM_JOB_ID}          # Job ID
msg 'Job name ' $SLURM_JOB_NAME        # Job name

################################################################################
module load R/4.5.1
module load pandoc/3.2
module load gdal/3.9.0
module load ImageMagick/7.1.1
module load gcc/14.2

Rscript -e "
rmarkdown::render(
  input = 'merscope-human-breast-cancer.Rmd',
  output_dir = '../docs',
  knit_root_dir = '.',
  output_format = rmarkdown::html_document(
    toc = TRUE,
    toc_float = TRUE,
    self_contained = FALSE,
    lib_dir = '../docs/site_libs',
    theme = 'cosmo'
  )
)
"


echo "------------------------------------------------------------------------"
echo "Job Completed on $(date)"
echo "------------------------------------------------------------------------"
