#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Utils modules
include { extractFastqPairFromDir } from './NextflowModules/Utils/fastq.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

// Mapping modules
include { BWAMapping } from './NextflowModules/BWA-Mapping/bwa-0.7.17_samtools-1.9/Mapping.nf' params(
    genome_fasta: "$params.genome", optional: '-c 100 -M'
)
include { MarkdupMerge  as Sambamba_MarkdupMerge } from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// IndelRealignment modules
include { RealignerTargetCreator as GATK_RealignerTargetCreator } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/RealignerTargetCreator.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_rtc_options"
)
include { IndelRealigner as GATK_IndelRealigner } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_ir_options"
)
include { ViewUnmapped as Sambamba_ViewUnmapped } from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include { Merge as Sambamba_Merge } from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// HaplotypeCaller modules
include { IntervalListTools as PICARD_IntervalListTools } from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count: "500", optional: ""
)
include { HaplotypeCaller as GATK_HaplotypeCaller } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_hc_options"
)
include { VariantFiltrationSnpIndel as GATK_VariantFiltration } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", snp_filter: "$params.gatk_snp_filter",
    snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)
include { CombineVariants as GATK_CombineVariants } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "--assumeIdenticalSamples"
)
include { SelectVariantsSample as GATK_SingleSampleVCF } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/SelectVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome"
)

// Fingerprint modules
include { UnifiedGenotyper as GATK_UnifiedGenotyper } from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome",
    optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES"
)

// QC Modules
include { FastQC } from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include { CollectMultipleMetrics as PICARD_CollectMultipleMetrics } from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(
    genome: "$params.genome",
    optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include { EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity } from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(
    optional: "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
)
include { CollectHsMetrics as PICARD_CollectHsMetrics } from './NextflowModules/Picard/2.22.0/CollectHsMetrics.nf' params(
    genome: "$params.genome", bait:"$params.dxtracks_path/$params.picard_bait",
    target: "$params.dxtracks_path/$params.picard_target",
    optional: "METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include { Flagstat as Sambamba_Flagstat } from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include { MultiQC } from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(
    optional: "--config $projectDir/assets/multiqc_config.yaml"
)
include { Mosdepth } from './NextflowModules/Mosdepth/0.3.3/Mosdepth.nf' params(optional: "-n -b $params.dxtracks_path/$params.exoncov_bed -Q 20")
include { Subsample as Sambamba_ViewSubsample } from './NextflowModules/Sambamba/0.7.0/ViewSubsample.nf' params(optional: "-L ${params.contamination_sites_bed}")
include { VerifyBamID2 } from './NextflowModules/VerifyBamID/2.0.1--h32f71e1_2/VerifyBamID2.nf'

//SNParray-calling Modules
include { IntervalListTools as PICARD_IntervalListToolsSNP } from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count:"100", optional: ""
)
include { HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF } from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(
    genome:"$params.genome", emit_ref_confidence: "BP_RESOLUTION", compress:true, optional: ""
)
include { GenotypeGVCF as GATK_GenotypeGVCF } from "./NextflowModules/GATK/4.2.1.0/GenotypeGvcfs.nf" params(
    genome:"$params.genome",  optional: "-all-sites", compress:true
)
include { MergeVcfs as GATK_MergeVcfs } from './NextflowModules/GATK/4.2.1.0/MergeVcfs.nf' params(
    genome:"$params.genome", compress:true
)

// CustomModules
include { IGV as BAF_IGV } from './CustomModules/BAF/IGV.nf'
include { CheckQC } from './CustomModules/CheckQC/CheckQC.nf'
include { SampleIndications as ClarityEpp_SampleIndications } from './CustomModules/ClarityEpp/SampleIndications.nf'
include { CallCNV as ExomeDepth_CallCNV } from './CustomModules/ExomeDepth/CallCNV.nf'
include { GetRefset as ExomeDepth_GetRefset } from './CustomModules/ExomeDepth/GetRefset.nf'
include { SingleIGV as ExomeDepth_SingleIGV } from './CustomModules/ExomeDepth/IGV.nf'
include { FamilyIGV as ExomeDepth_FamilyIGV } from './CustomModules/ExomeDepth/IGV.nf'
include { Summary as ExomeDepth_Summary } from './CustomModules/ExomeDepth/Summary.nf'
include { ImportBam as ExonCov_ImportBam } from './CustomModules/ExonCov/ImportBam.nf'
include { SampleQC as ExonCov_SampleQC } from './CustomModules/ExonCov/SampleQC.nf'
include { IGV as UPD_IGV } from './CustomModules/UPD/IGV.nf'
include { CreateHSmetricsSummary } from './CustomModules/Utils/CreateHSmetricsSummary.nf'
include { GetStatsFromFlagstat } from './CustomModules/Utils/GetStatsFromFlagstat.nf'
include { Kinship } from './CustomModules/Kinship/Kinship.nf'
include { ParseChildFromFullTrio } from './CustomModules/Utils/ParseChildFromFullTrio.nf'
include { SavePedFile } from './CustomModules/Utils/SavePedFile.nf'
include { VersionLog } from './CustomModules/Utils/VersionLog.nf'
include { Fraction } from './CustomModules/Utils/ParseDownsampleFraction.nf'


def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter GATK_RealignerTargetCreator
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

// Define ped file, used in Kinship
def ped_file = file("${params.ped_folder}/${analysis_id}.ped")
if (!ped_file.exists()) {
    error("ERROR: ${ped_file} not found.")
}

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

    // GATK HaplotypeCaller
    PICARD_IntervalListTools(Channel.fromPath("$params.dxtracks_path/$params.gatk_hc_interval_list"))
    GATK_HaplotypeCaller(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, bam_file, bai_file]}
            .groupTuple()
            .combine(PICARD_IntervalListTools.out.flatten())
    )
    GATK_VariantFiltration(GATK_HaplotypeCaller.out)
    GATK_CombineVariants(GATK_VariantFiltration.out.groupTuple())
    GATK_SingleSampleVCF(GATK_CombineVariants.out.combine(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id]})
    )

    // GATK HaplotypeCaller (SNParray target)
    PICARD_IntervalListToolsSNP(Channel.fromPath("$params.dxtracks_path/$params.gatk_hc_interval_list_snparray"))
    GATK_HaplotypeCallerGVCF(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}
        .groupTuple()
        .combine(PICARD_IntervalListToolsSNP.out.flatten())
    )
    GATK_GenotypeGVCF(GATK_HaplotypeCallerGVCF.out)
    GATK_MergeVcfs(GATK_GenotypeGVCF.out.groupTuple())

    //Downsampling BAM files for VerifyBAMID
    Mosdepth(Sambamba_Merge.out)
    Fraction(Mosdepth.out.summary_file)
    Sambamba_ViewSubsample(Fraction.out.join(Sambamba_Merge.out))

    // ExomeDepth
    ExomeDepth_CallCNV(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})
    ExomeDepth_Summary(analysis_id, ExomeDepth_CallCNV.out.HC_stats_log.collect())

    // ExomeDepth IGV single sessions
    ExomeDepth_GetRefset(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file]}.groupTuple())
    ExomeDepth_SingleIGV(ExomeDepth_GetRefset.out.map{sample_id, refset -> [sample_id, analysis_id, refset]})

    // BAF analysis per sample
    BAF_IGV(GATK_MergeVcfs.out)

    // UPD analysis per family
    ParseChildFromFullTrio(ped_file, GATK_MergeVcfs.out.map{sample_id, vcf_file, vcf_idx_file -> [sample_id]}.collect())
    UPD_IGV(
        ped_file,
        analysis_id,
        ParseChildFromFullTrio.out.splitCsv().flatten(),
        GATK_MergeVcfs.out.map{output_name, vcf_files, vcf_idx_files -> [vcf_files]}.collect()
    )

    //ExomeDepth IGV family sessions
    ExomeDepth_FamilyIGV(
        ped_file,
        analysis_id,
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file]}.collect(),
        GATK_SingleSampleVCF.out.map{sample_id, vcf_files, vcf_idx_files -> [vcf_files]}.collect(),
        ExomeDepth_CallCNV.out.HC_vcf.collect(),
        ExomeDepth_CallCNV.out.HC_igv.collect(),
        UPD_IGV.out.UPD_IGV_files.collect(),
        BAF_IGV.out.BAF_IGV_files.collect()
    )

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper(Sambamba_Merge.out)

    // QC - FastQC
    FastQC(fastq_files)

    // QC - Picard
    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_Merge.out)
    PICARD_CollectHsMetrics(Sambamba_Merge.out)
    CreateHSmetricsSummary(PICARD_CollectHsMetrics.out.collect())

    // QC - Flagstat
    Sambamba_Flagstat(Sambamba_Merge.out)
    GetStatsFromFlagstat(Sambamba_Flagstat.out.collect())

    // QC - ExonCov
    ExonCov_ImportBam(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]}
    )
    ClarityEpp_SampleIndications(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> sample_id})
    ExonCov_SampleQC(
        ExonCov_ImportBam.out.join(ClarityEpp_SampleIndications.out)
            .map{sample_id, exoncov_id, indication -> [analysis_id, exoncov_id, indication]}
            .groupTuple()
    )

    // QC - VerifyBamID2 (contamination)
    VerifyBamID2(Sambamba_ViewSubsample.out.groupTuple())

    // QC - MultiQC report
    MultiQC(analysis_id, Channel.empty().mix(
        FastQC.out,
        PICARD_CollectMultipleMetrics.out,
        PICARD_EstimateLibraryComplexity.out,
        PICARD_CollectHsMetrics.out,
        VerifyBamID2.out.map{sample_id, self_sm -> [self_sm]},
        ExonCov_SampleQC.out
    ).collect())

    // QC - Kinship
    Kinship(GATK_CombineVariants.out, ped_file)

    // QC - Check and collect
    CheckQC(
        analysis_id,
        Channel.empty().mix(
            MultiQC.out.map{html, report_data_dir -> [report_data_dir + '/multiqc_picard_HsMetrics.txt']},
            MultiQC.out.map{html, report_data_dir -> [report_data_dir + '/multiqc_verifybamid.txt']},
            Kinship.out.map{analysis, kinship_name, kinship_check_out -> [kinship_check_out]},
        ).collect()
    )

    //SavePedFile
    SavePedFile(ped_file)

    // Create log files: Repository versions and Workflow params
    VersionLog(Channel.of(
        "${workflow.projectDir}/",
        "${params.dxtracks_path}/",
        "${params.exoncov_path}/",
        "${params.clarity_epp_path}/",
        "${params.dx_resources_path}/",
        "${params.upd_path}/",
        "${params.baf_path}/",
    ).collect())
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$projectDir/assets/workflow_complete.html")
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
