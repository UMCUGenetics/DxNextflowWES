#!/bin/bash
# Set input and output dirs
input=`realpath $1`
output=`realpath $2`
email=$3
mkdir -p $output && cd $output

workflow_path='/hpc/diaggen/software/development/DxNextflowWES'

sbatch <<EOT
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH -o slurm_nextflow_wes.%j.out
#SBATCH -e slurm_nextflow_wes.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/WES.nf \
-c $workflow_path/WES.config \
--fastq_path $input \
--outdir $output \
--email $email \
-profile slurm \
-resume -ansi-log false
EOT
