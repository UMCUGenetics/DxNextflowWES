#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "$params.gatk_known_indels")

include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'

fastq_files = extractFastqPairFromDir(params.fastq_path)

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

    // ToDo:
    // Mapping
        // bwa mem
        // sambamba sort, merge + markdup

    // QC
        // FastQC
        // flagstat
        // bammetrics replacement

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
