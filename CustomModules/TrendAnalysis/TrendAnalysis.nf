process TrendAnalysis {
    // Custom process to run Trend_Analysis_tool
    tag {"TrendAnalysis ${analysis_id}"}
    label 'TrendAnalysis'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, path(input_files))

    script:
        """
        source ${params.trend_analysis_path}/venv/bin/activate
        python ${params.trend_analysis_path}/trend_analysis.py upload processed_data ${analysis_id} .
        """
}