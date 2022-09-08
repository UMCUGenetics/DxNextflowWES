#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

// Mapping modules
include BWAMapping from './NextflowModules/BWA-Mapping/bwa-0.7.17_samtools-1.9/Mapping.nf' params(
    genome_fasta: "$params.genome", optional: '-c 100 -M'
)
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// IndelRealignment modules
include RealignerTargetCreator as GATK_RealignerTargetCreator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/RealignerTargetCreator.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_rtc_options"
)
include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: ""
)
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome",
    optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES"
)

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(
    genome: "$params.genome",
    optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(
    optional: "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
)
include CollectHsMetrics as PICARD_CollectHsMetrics from './NextflowModules/Picard/2.22.0/CollectHsMetrics.nf' params(
    genome: "$params.genome", bait:"$params.dxtracks_path/$params.picard_bait",
    target: "$params.dxtracks_path/$params.picard_target",
    optional: "METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include MultiQC from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(
    optional: "--config $baseDir/assets/multiqc_config.yaml"
)
include VerifyBamID2 from './NextflowModules/VerifyBamID/2.0.1--h32f71e1_2/VerifyBamID2.nf'

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter GATK_RealignerTargetCreator
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

workflow {
    // Mapping
    BWAMapping(fastq_files)
    Sambamba_MarkdupMerge(
        BWAMapping.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file]
        }.groupTuple()
    )

    // GATK IndelRealigner
    GATK_RealignerTargetCreator(Sambamba_MarkdupMerge.out.combine(chromosomes))
    GATK_IndelRealigner(Sambamba_MarkdupMerge.out.combine(GATK_RealignerTargetCreator.out, by: 0))
    Sambamba_ViewUnmapped(Sambamba_MarkdupMerge.out)
    Sambamba_Merge(GATK_IndelRealigner.out.mix(Sambamba_ViewUnmapped.out).groupTuple())

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper(Sambamba_Merge.out)

    // Clarity epp
    ClarityEppIndications(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> sample_id})

    // ExonCov
    ExonCovImportBam(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]}
    )

    // QC
    FastQC(fastq_files)

    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_Merge.out)
    PICARD_CollectHsMetrics(Sambamba_Merge.out)
    CreateHSmetricsSummary(PICARD_CollectHsMetrics.out.collect())

    Sambamba_Flagstat(Sambamba_Merge.out)
    GetStatsFromFlagstat(Sambamba_Flagstat.out.collect())

    ExonCovSampleQC(
        ExonCovImportBam.out.join(ClarityEppIndications.out)
            .map{sample_id, exoncov_id, indication -> [analysis_id, exoncov_id, indication]}
            .groupTuple()
    )

    VerifyBamID2(Sambamba_Merge.out.groupTuple())

    MultiQC(analysis_id, Channel.empty().mix(
        FastQC.out,
        PICARD_CollectMultipleMetrics.out,
        PICARD_EstimateLibraryComplexity.out,
        PICARD_CollectHsMetrics.out,
        VerifyBamID2.out.map{sample_id, self_sm -> [self_sm]},
        ExonCovSampleQC.out
    ).collect())

    // Create log files: Repository versions and Workflow params
    VersionLog()
    Workflow_ExportParams()
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
        def subject = "WES Preprocess Workflow Successful: ${analysis_id}"
        sendMail(
            to: params.email.trim(),
            subject: subject,
            body: email_html,
            attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html"
        )

    } else {
        def subject = "WES Preprocess Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}

// Custom processes
process ExonCovImportBam {
    // Custom process to run ExonCov import_bam
    tag {"ExonCov ImportBam ${sample_id}"}
    label 'ExonCov'
    label 'ExonCov_ImportBam'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py import_bam \
        --threads ${task.cpus} \
        --overwrite \
        --print_sample_id \
        --exon_bed ${params.dxtracks_path}/${params.exoncov_bed} \
        ${analysis_id} WES ${bam_file} | tr -d '\n'
        """
}

process ExonCovSampleQC {
    // Custom process to run ExonCov sample_qc
    tag {"ExonCov Sample QC ${analysis_id}"}
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

process ClarityEppIndications {
    // Custom process to run clarity_epp export sample_indications
    tag {"ClarityEppExportSampleIndications ${sample_id}"}
    label 'ClarityEpp'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a clarity export restarting the workflow.

    input:
        val(sample_id)

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.clarity_epp_path}/venv/bin/activate
        python ${params.clarity_epp_path}/clarity_epp.py export sample_indications \
        -a ${sample_id} | cut -f 2 | grep -v 'Indication' | tr -d '\n'
        """
}

process GetStatsFromFlagstat {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"GetStatsFromFlagstat"}
    label 'GetStatsFromFlagstat'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(flagstat_files)

    output:
        path('run_stats.txt')

    script:
        """
        python ${baseDir}/assets/get_stats_from_flagstat.py ${flagstat_files} > run_stats.txt
        """
}

process CreateHSmetricsSummary {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"CreateHSmetricsSummary"}
    label 'CreateHSmetricsSummary'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hsmetrics_files)

    output:
        path('HSMetrics_summary.txt')

    script:
        """
        python ${baseDir}/assets/create_hsmetrics_summary.py ${hsmetrics_files} > HSMetrics_summary.txt
        """
}

process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowWes' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_tracks' >> repository_version.log
        git --git-dir=${params.dxtracks_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'ExonCov' >> repository_version.log
        git --git-dir=${params.exoncov_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'clarity_epp' >> repository_version.log
        git --git-dir=${params.clarity_epp_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
