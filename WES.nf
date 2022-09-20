#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include extractBamFromDir from './NextflowModules/Utils/bam.nf'
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

// HaplotypeCaller modules
include IntervalListTools as PICARD_IntervalListTools from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count: "500", optional: ""
)
include HaplotypeCaller as GATK_HaplotypeCaller from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_hc_options"
)
include VariantFiltrationSnpIndel as GATK_VariantFiltration from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", snp_filter: "$params.gatk_snp_filter",
    snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)
include CombineVariants as GATK_CombineVariants from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--assumeIdenticalSamples"
)
include SelectVariantsSample as GATK_SingleSampleVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/SelectVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome"
)

//SNParray-calling Modules
include IntervalListTools as PICARD_IntervalListToolsSNP from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count:"100", optional: ""
)
include HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(
    genome:"$params.genome", emit_ref_confidence: "BP_RESOLUTION", compress:true, optional: ""
)
include GenotypeGVCF as GATK_GenotypeGVCF from "./NextflowModules/GATK/4.2.1.0/GenotypeGvcfs.nf" params(
    genome:"$params.genome",  optional: "-all-sites", compress:true
)
include MergeVcfs as GATK_MergeVcfs from './NextflowModules/GATK/4.2.1.0/MergeVcfs.nf' params(
    genome:"$params.genome", compress:true
)

// CustomModules
include IGV as BAF_IGV from './CustomModules/BAF/IGV.nf'
include CallCNV as ExomeDepth_CallCNV from './CustomModules/ExomeDepth/CallCNV.nf'
include GetRefset as ExomeDepth_GetRefset from './CustomModules/ExomeDepth/GetRefset.nf'
include SingleIGV as ExomeDepth_SingleIGV from './CustomModules/ExomeDepth/IGV.nf'
include FamilyIGV as ExomeDepth_FamilyIGV from './CustomModules/ExomeDepth/IGV.nf'
include Summary as ExomeDepth_Summary from './CustomModules/ExomeDepth/ExomeDepthSummary.nf'
include IGV as UPD_IGV from './CustomModules/UPD/IGV.nf'
include Kinship from './CustomModules/Utils/Kinship.nf'
include ParseChildFromFullTrio from './CustomModules/Utils/ParseChildFromFullTrio.nf'
include VersionLog from './CustomModules/Utils/VersionLog.nf'

def input_bam_files = extractBamFromDir(params.bam_path)
def analysis_id = params.outdir.split('/')[-1]

// Define ped file, used in Kinship
def ped_file = file("${params.ped_folder}/${analysis_id}.ped")
if (!ped_file.exists()) {
    exit 1, "ERROR: ${ped_file} not found."
}

workflow {
    // GATK HaplotypeCaller
    PICARD_IntervalListTools(Channel.fromPath("$params.dxtracks_path/$params.gatk_hc_interval_list"))
    GATK_HaplotypeCaller(
        input_bam_files.map{sample_id, bam_file, bai_file -> [analysis_id, bam_file, bai_file]}
            .groupTuple()
            .combine(PICARD_IntervalListTools.out.flatten())
    )
    GATK_VariantFiltration(GATK_HaplotypeCaller.out)
    GATK_CombineVariants(GATK_VariantFiltration.out.groupTuple())
    GATK_SingleSampleVCF(GATK_CombineVariants.out.combine(
        input_bam_files.map{sample_id, bam_file, bai_file -> [sample_id]})
    )

    // ExomeDepth
    ExomeDepth_CallCNV(input_bam_files.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})
    ExomeDepth_Summary(analysis_id, ExomeDepth_CallCNV.out.HC_stats_log.collect())

    // ExomeDepth IGV sessions
    ExomeDepth_GetRefset(input_bam_files.map{sample_id, bam_file, bai_file -> [sample_id, bam_file]}.groupTuple())
    ExomeDepth_SingleIGV(ExomeDepth_GetRefset.out.map{sample_id, refset -> [sample_id, analysis_id, refset]})
    ParseChildFromFullTrio(ped_file, GATK_MergeVcfs.out.map{sample_id, vcf_file, vcf_idx_file -> [sample_id]}.collect())
    ExomeDepth_FamilyIGV(ParseChildFromFullTrio.out.splitCsv().flatten()
        .combine(input_bam_files.map{ sample_id, bam_file, bai_file -> [bam_file]}).groupTuple()
        .map{sample_id, bam_files -> [sample_id, bam_files, ped_file, analysis_id]}
    )

    // Kinship
    Kinship(GATK_CombineVariants.out)

    // GATK HaplotypeCaller (SNParray target)
    PICARD_IntervalListToolsSNP(Channel.fromPath("$params.dxtracks_path/$params.gatk_hc_interval_list_snparray"))
    GATK_HaplotypeCallerGVCF(
        input_bam_files.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}
        .groupTuple()
        .combine(PICARD_IntervalListToolsSNP.out.flatten())
    )
    GATK_GenotypeGVCF(GATK_HaplotypeCallerGVCF.out)
    GATK_MergeVcfs(GATK_GenotypeGVCF.out.groupTuple())

    // BAF analysis per sample
    BAF_IGV(GATK_MergeVcfs.out)

    // UPD analysis per family
    UPD_IGV(
        ped_file,
        analysis_id,
        ParseChildFromFullTrio.out.splitCsv().flatten(),
        GATK_MergeVcfs.out.map{output_name, vcf_files, vcf_idx_files -> [vcf_files]}.collect()
    )

    //SavePedFile
    SavePedFile()

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
        def subject = "WES Workflow Successful: ${analysis_id}"
        sendMail(
            to: params.email.trim(),
            subject: subject,
            body: email_html,
            attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html"
        )

    } else {
        def subject = "WES Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
