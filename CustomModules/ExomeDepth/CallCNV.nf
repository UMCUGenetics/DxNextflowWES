process CallCNV {
    // Custom process to run Exomedepth
    tag {"ExomeDepth CallCNV ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_CallCNV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    output:
        path("*.log", emit: ED_log)
        path("HC_*.igv", emit: HC_igv)
        path("UMCU_*.igv", emit: UMCU_igv)
        path("HC_*.vcf", emit: HC_vcf)
        path("UMCU_*.vcf", emit: UMCU_vcf)
        path('HC_*_stats.log', emit: HC_stats_log)
        path('UMCU_*_stats.log', emit: UMCU_stats_log)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/run_ExomeDepth.py callcnv ./ ${bam_file} ${analysis_id} ${sample_id}
        """
}