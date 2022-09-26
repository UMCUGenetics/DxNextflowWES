process GetRefset{
    // Custom process to run get exomedepth reference set from exomedepth db
    tag {"ExomeDepth GetRefset ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_GetRefset'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false

    input:
        tuple(sample_id, path(bam_file))

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/exomedepth_db.py add_sample_return_refset_bam ${bam_file} --print_refset_stdout | tr -d '\n'
        """
}