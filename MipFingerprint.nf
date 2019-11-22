#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf'
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(params)
include MipsTrimDedup from './NextflowModules/Mips/1.0.1/MipsTrimDedup.nf' params(
    outdir: params.outdir,
    mips_trim_dedup: params.mips_trim_dedup,
    design_file: params.mip_design_file,
    uuid_length: params.mip_uuid_length,
    uuid_read: params.mip_uuid_read,
    )
include MEM as BWAMEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(params)
include ViewSort as SambambaViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf' params(params)
include UnifiedGenotyper as GATKUnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(params)

fastq_files = extractFastqFromDir(params.fastq_path)
samples = fastq_files.groupTuple(by:[0])

workflow MipFingerprint {
    FastQC(fastq_files)
    MipsTrimDedup(samples)
    BWAMEM(MipsTrimDedup.out)
    SambambaViewSort(BWAMEM.out)
    GATKUnifiedGenotyper(SambambaViewSort.out)
}
