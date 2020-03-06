#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MipsTrimDedup from './NextflowModules/Mips/1.0.1/MipsTrimDedup.nf' params(trim_dedup_path: "$params.mips_trim_dedup_path", design_file: "$params.mips_design_file", uuid_length: "$params.mips_uuid_length", uuid_read: "$params.mips_uuid_read")
include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--intervals $params.fingerprint_target --output_mode EMIT_ALL_SITES")
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:'')

fastq_files = extractFastqPairFromDir(params.fastq_path)
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

    // ToDo:
    // cleanup script -> QC, extra log files?
}
