#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_ir_options")
include HaplotypeCaller as GATK_HaplotypeCaller from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_hc_options")

include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics")
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf'
include CollectHsMetrics as PICARD_CollectHsMetrics from './NextflowModules/Picard/2.22.0/CollectHsMetrics.nf' params(genome:"$params.genome", bait:"$params.picard_bait", target:"$params.picard_target", optional: "METRIC_ACCUMULATION_LEVEL=SAMPLE")
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:'')

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Mapping
    BWA_MEM(fastq_files)
    Sambamba_ViewSort(BWA_MEM.out)
    Sambamba_MarkdupMerge(
        Sambamba_ViewSort.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file]
        }.groupTuple()
    )

    // GATk
    GATK_IndelRealigner(Sambamba_MarkdupMerge.out)


    // QC
    FastQC(fastq_files)
    Sambamba_Flagstat(Sambamba_MarkdupMerge.out)
    PICARD_CollectMultipleMetrics(Sambamba_MarkdupMerge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_MarkdupMerge.out)
    PICARD_CollectHsMetrics(Sambamba_MarkdupMerge.out)

    MultiQC(Channel.empty().mix(
        FastQC.out,
        Sambamba_Flagstat.out,
        PICARD_CollectMultipleMetrics.out,
        PICARD_EstimateLibraryComplexity.out,
        PICARD_CollectHsMetrics.out
    ).collect())

    // ToDo:
    // Mapping
        // bwa mem
        // sambamba sort, merge + markdup

    // QC
        // FastQC
        // flagstat
        // bammetrics replacement
            // Picard CollectMultipleMetrics
                // PROGRAM=CollectAlignmentSummaryMetrics -> multiqc
                // PROGRAM=CollectInsertSizeMetrics -> multiqc
                //PROGRAM=QualityScoreDistribution
            // Picard EstimateLibraryComplexity
            // Picard CalculateHsMetrics -> multiqc

    // GATk
        // realignment
        // haplotypecaller
        // filter
        // fingerprint

    // ExonCov

    // Other
        // Kinship
        // Single sample vcf
        // Gendercheck
        // cleanup

}
