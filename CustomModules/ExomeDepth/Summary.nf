process Summary {
    // Custom process to stats from ExomeDepth analysis
    tag {"ExomeDepth Summary ${analysis_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_Summary'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(analysis_id)
        path(exomedepth_logs)

    output:
        path("${analysis_id}_exomedepth_summary.txt")

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/run_ExomeDepth.py summary ${exomedepth_logs} > ${analysis_id}_exomedepth_summary.txt
        """
}