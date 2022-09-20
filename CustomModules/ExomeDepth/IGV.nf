process SingleIGV {
    // Custom process to run Single sample IGV analysis
    tag {"ExomeDepth SingleIGV ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_SingleIGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, analysis_id, refset)

    output:
        path("*.xml", emit: Single_IGV_file)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/igv_xml_session.py single_igv ./ ${sample_id} ${analysis_id} ${refset}
        """
}

process FamilyIGV {
    // Custom process to run Family IGV analysis
    tag {"ExomeDepth FamilyIGV ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_FamilyIGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, path(bam_files), ped_file, analysis_id)

    output:
        path("*.xml", emit: Family_IGV_file)

    script:
        def bam_files = bam_files.collect{"$it"}.join(" ")
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/igv_xml_session.py family_igv ./ ${ped_file} ${analysis_id} ${sample_id} ${bam_files}
        """
}