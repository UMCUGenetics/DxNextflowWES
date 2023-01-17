process SampleQC {
    // Custom process to run ExonCov sample_qc
    tag {"ExonCov SampleQC ${analysis_id}"}
    label 'ExonCov'
    label 'ExonCov_SampleQC'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_ids, indications)

    output:
        path("${analysis_id}.ExonCovQC_check.out")

    script:
        def samples = sample_ids.collect{"$it"}.join(" ")
        def panels = indications.collect{"$it"}.join(" ")
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py sample_qc \
        -s ${samples} -p ${panels} > ${analysis_id}.ExonCovQC_check.out
        """
}