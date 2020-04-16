#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// Mapping modules
include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// IndelRealignment modules
include RealignerTargetCreator as GATK_RealignerTargetCreator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/RealignerTargetCreator.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_rtc_options")
include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "")
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// HaplotypeCaller modules
include IntervalListTools as PICARD_IntervalListTools from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(interval_list: "$params.gatk_hc_interval_list", scatter_count:'500')
include HaplotypeCaller as GATK_HaplotypeCaller from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_hc_options")
include VariantFiltrationSnpIndel as GATK_VariantFiltration from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome:"$params.genome", snp_filter: "$params.gatk_snp_filter", snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)
include CombineVariants as GATK_CombineVariants from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--assumeIdenticalSamples")
include SelectVariantsSample as GATK_SingleSampleVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/SelectVariants.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome")

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--intervals $params.fingerprint_target --output_mode EMIT_ALL_SITES")

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(optional:"OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500")
include CollectHsMetrics as PICARD_CollectHsMetrics from './NextflowModules/Picard/2.22.0/CollectHsMetrics.nf' params(genome:"$params.genome", bait:"$params.picard_bait", target:"$params.picard_target", optional: "METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:"--config $baseDir/assets/multiqc_config.yaml")

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter GATK_RealignerTargetCreator
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

// Define ped file, used in Kinship
def ped_file = file("${params.ped_folder}/${analysis_id}.ped")
if (!ped_file.exists()) {
    exit 1, "${ped_file.getName} not found in ${ped_file.getParent}."
}

workflow {
    // Mapping
    BWA_MEM(fastq_files)
    Sambamba_ViewSort(BWA_MEM.out)
    Sambamba_MarkdupMerge(
        Sambamba_ViewSort.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file]
        }.groupTuple()
    )

    // GATK IndelRealigner
    GATK_RealignerTargetCreator(Sambamba_MarkdupMerge.out.combine(chromosomes))
    GATK_IndelRealigner(Sambamba_MarkdupMerge.out.combine(GATK_RealignerTargetCreator.out, by: 0))
    Sambamba_ViewUnmapped(Sambamba_MarkdupMerge.out)
    Sambamba_Merge(GATK_IndelRealigner.out.mix(Sambamba_ViewUnmapped.out).groupTuple())

    // GATK HaplotypeCaller
    PICARD_IntervalListTools(Channel.fromPath(params.gatk_hc_interval_list))
    GATK_HaplotypeCaller(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, bam_file, bai_file]}.groupTuple().combine(PICARD_IntervalListTools.out.flatten()))
    GATK_VariantFiltration(GATK_HaplotypeCaller.out)
    GATK_CombineVariants(GATK_VariantFiltration.out.groupTuple())
    GATK_SingleSampleVCF(GATK_CombineVariants.out.combine(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id]}))

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper(Sambamba_Merge.out)

    // ExonCov
    ExonCov(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})

    // Kinship
    Kinship(GATK_CombineVariants.out)

    // QC
    FastQC(fastq_files)

    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_Merge.out)
    PICARD_CollectHsMetrics(Sambamba_Merge.out)
    CreateHSmetricsSummary(PICARD_CollectHsMetrics.out.collect())

    Sambamba_Flagstat(Sambamba_Merge.out)
    GetStatsFromFlagstat(Sambamba_Flagstat.out.collect())

    MultiQC(
        Channel.empty().mix(
            FastQC.out.flatten().map{file -> [analysis_id, file]},
            PICARD_CollectMultipleMetrics.out.flatten().map{file -> [analysis_id, file]},
            PICARD_EstimateLibraryComplexity.out.map{file -> [analysis_id, file]},
            PICARD_CollectHsMetrics.out.map{file -> [analysis_id, file]}
        ).groupTuple()
    )

    TrendAnalysisTool(
        GATK_CombineVariants.out.map{id, vcf_file, idx_file -> [id, vcf_file]}
            .concat(GetStatsFromFlagstat.out.map{file -> [analysis_id, file]})
            .concat(CreateHSmetricsSummary.out.map{file -> [analysis_id, file]})
            .groupTuple()
    )
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "WES Workflow Successful: ${analysis_id}"
        sendMail(to: params.email, subject: subject, body: email_html, attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "WES Workflow Failed: ${analysis_id}"
        sendMail(to: params.email, subject: subject, body: email_html)
    }
}

// Custom processes
process ExonCov {
    // Custom process to run ExonCov
    tag {"ExonCov ${sample_id}"}
    label 'ExonCov'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
    tuple analysis_id, sample_id, file(bam_file), file(bai_file)

    script:
    """
    source ${params.exoncov_path}/venv/bin/activate
    python ${params.exoncov_path}/ExonCov.py import_bam --threads ${task.cpus} --overwrite --exon_bed ${params.exoncov_bed} ${analysis_id} ${bam_file}
    """
}

process Kinship {
    // Custom process to run Kinship tools
    // Container does not work
    // king: run.c:355: main: Unexpected error: No such file or directory.
    // Aborted (core dumped)
    tag {"Kinship ${analysis_id}"}
    label 'Kinship'
    //container = '/hpc/diaggen/software/guix_containers/kinship.sif'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    tuple analysis_id, file(vcf_file), file(vcf_index)

    output:
    tuple analysis_id, file("${analysis_id}.kinship"), file("${analysis_id}.kinship_check.out")

    script:
    """
    ${params.vcftools_path}/vcftools --vcf ${vcf_file} --plink
    ${params.plink_path}/plink --file out --make-bed --noweb
    ${params.king_path}/king -b plink.bed --kinship
    cp king.kin0 ${analysis_id}.kinship
    python ${baseDir}/assets/check_kinship.py ${analysis_id}.kinship ${ped_file} > ${analysis_id}.kinship_check.out
    """
}

process GetStatsFromFlagstat {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"GetStatsFromFlagstat"}
    label 'GetStatsFromFlagstat'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
    file(flagstat_files: "*")

    output:
    file('run_stats.txt')

    script:
    """
    python ${baseDir}/assets/get_stats_from_flagstat.py ${flagstat_files} > run_stats.txt
    """
}

process CreateHSmetricsSummary {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"CreateHSmetricsSummary"}
    label 'CreateHSmetricsSummary'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
    file(hsmetrics_files: "*")

    output:
    file('HSMetrics_summary.txt')

    script:
    """
    python ${baseDir}/assets/create_hsmetrics_summary.py ${hsmetrics_files} > HSMetrics_summary.txt
    """
}

process TrendAnalysisTool {
    // Custom process to run Trend_Analysis_tool
    tag {"TrendAnalysisTool ${analysis_id}"}
    label 'TrendAnalysisTool'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
    tuple analysis_id, file(input_files: "*")

    script:
    """
    source ${params.trend_analysis_path}/venv/bin/activate
    python ${params.trend_analysis_path}/trend_analysis.py upload processed_data ${analysis_id} .
    """
}
