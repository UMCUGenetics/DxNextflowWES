process ParseChildFromFullTrio {
    //Custom process to parse PED file and output sampleID of children with both parents.
    tag {"ParseChildFromFullTrio ${analysis_id}"}
    label 'ParseChildFromFullTrio'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false

    input:
        path(ped_file)
        val(sample_id)

    output:
        stdout emit: trio_sample

    script:
        def sample_ids = sample_id.join(" ")
        """
        python ${baseDir}/assets/parse_child_from_fulltrio.py ${ped_file} ${sample_ids} | tr -d '\n'
        """
}