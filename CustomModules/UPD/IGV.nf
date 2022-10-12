process IGV {
    // Custom process to run UPD analysis
    tag {"UPD IGV $trio_sample"}
    label 'UPD'
    label 'UPD_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(ped_file)
        val(analysis_id)
        val(trio_sample)
        path(vcf_files)

    output:
        path("*.igv", emit: UPD_IGV_files)

    script:
        """
        source ${params.upd_path}/venv/bin/activate
        python ${params.upd_path}/make_UPD_igv.py ${ped_file} ${analysis_id} $trio_sample ${vcf_files} -c
        """
}