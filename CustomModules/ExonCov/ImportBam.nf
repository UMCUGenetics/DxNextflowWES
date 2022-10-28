process ImportBam {
    // Custom process to run ExonCov import_bam
    tag {"ExonCov ImportBam ${sample_id}"}
    label 'ExonCov'
    label 'ExonCov_ImportBam'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py import_bam \
        --threads ${task.cpus} \
        --overwrite \
        --print_sample_id \
        --exon_bed ${params.dxtracks_path}/${params.exoncov_bed} \
        ${analysis_id} WES ${bam_file} | tr -d '\n'
        """
}