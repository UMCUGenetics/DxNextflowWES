#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/production/DxNextflowWES'

# Set input and output dirs
input=`realpath -e $1`
output=`realpath $2`
email=$3

design=${4-default}
if [ $design == "SSv7" ]; then
    echo " #### SSv7 design target files are used for the analysis ####"
    gatk_hc_interval_list='--gatk_hc_interval_list Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v3_CREv2_SSv7_100bpflank.interval_list'
    picard_bait='--picard_bait Tracks/SureSelect_v7_elidS31285117_Covered.list'
fi

mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name Nextflow_WES
#SBATCH -o log/slurm_nextflow_wes.%j.out
#SBATCH -e log/slurm_nextflow_wes.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/WES.nf \
-c $workflow_path/WES.config \
--fastq_path $input \
--outdir $output \
--email $email \
${gatk_hc_interval_list:-""} \
${picard_bait:-""} \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    echo "Remove work directory"
    rm -r work

    echo "Creating md5sum"
    find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

    echo "WES workflow completed successfully."
    rm workflow.running
    touch workflow.done

    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi
