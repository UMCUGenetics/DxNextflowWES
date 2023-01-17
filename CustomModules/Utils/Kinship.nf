process Kinship {
    // Custom process to run Kinship tools
    tag {"Kinship ${analysis_id}"}
    label 'Kinship'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(analysis_id, path(vcf_file), path(vcf_index))
        path(ped_file)

    output:
        tuple(analysis_id, path("${analysis_id}.kinship"), path("${analysis_id}.kinship_check.out"))

    script:
        """
        ${params.vcftools_path}/vcftools --vcf ${vcf_file} --plink
        ${params.plink_path}/plink --file out --make-bed --noweb
        ${params.king_path}/king -b plink.bed --kinship
        cp king.kin0 ${analysis_id}.kinship
        python ${baseDir}/assets/check_kinship.py ${analysis_id}.kinship ${ped_file} > ${analysis_id}.kinship_check.out
        """
}