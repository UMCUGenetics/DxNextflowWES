#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/production/DxNextflowWES'

# Set input and output dirs
input=`realpath -e $1`
output=`realpath $2`
email=$3
mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name Nextflow_WES_Fingerprint
#SBATCH -o log/slurm_nextflow_wes_fingerprint.%j.out
#SBATCH -e log/slurm_nextflow_wes_fingerprint.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/WES_Fingerprint.nf \
-c $workflow_path/WES.config \
--bam_path $input \
--outdir $output \
--email $email \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    echo "Remove work directory"
    rm -r work

    echo "WES Fingerprint workflow completed successfully."
    rm workflow.running
    touch workflow.done

    echo "Change permissions"
    chmod 775 -R $output
    
    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    
    echo "Change permissions"
    chmod 775 -R $output
    
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi
