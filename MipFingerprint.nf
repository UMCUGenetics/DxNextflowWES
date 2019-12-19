#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf'
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(params)
include MipsTrimDedup from './NextflowModules/Mips/1.0.1/MipsTrimDedup.nf' params(params)
include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(params)
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf' params(params)
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf' params(params)
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(params)
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(params)

fastq_files = extractFastqFromDir(params.fastq_path)
samples = fastq_files.map( {it.flatten()}).groupTuple(by:[0])

workflow {
    FastQC(fastq_files)
    MipsTrimDedup(samples)
    BWA_MEM(MipsTrimDedup.out)
    Sambamba_ViewSort(BWA_MEM.out)
    Sambamba_Flagstat(Sambamba_ViewSort.out)
    GATK_UnifiedGenotyper(Sambamba_ViewSort.out)

    // Multi QC files
    multi_qc_files = Channel.empty().mix(FastQC.out, Sambamba_Flagstat.out).collect()
    MultiQC(multi_qc_files)
}
