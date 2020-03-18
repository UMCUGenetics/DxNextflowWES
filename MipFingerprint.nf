#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
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

// Custom processes
process MipsTrimDedup {
    // Custom process to run MIPS TrimDedup
    tag {"MIPS TrimDedup ${sample_id} - ${rg_id}"}
    label 'MIPS_1_0_1'
    label 'MIPS_1_0_1_TrimDedup'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    tuple sample_id, rg_id, file(r1_fastqs: "*"), file(r2_fastqs: "*")

    output:
    tuple sample_id, rg_id, file('*_LMergedTrimmedDedup_R1_*.fastq.gz'), file('*_LMergedTrimmedDedup_R2_*.fastq.gz')

    script:
    def r1_args = r1_fastqs.collect{ "$it" }.join(" ")
    def r2_args = r2_fastqs.collect{ "$it" }.join(" ")

    rg_id = "${sample_id}_MergedTrimmedDedup"

    """
    python ${params.trim_dedup_path} -d ${params.design_file}  -l ${params.uuid_length} -ur ${params.uuid_read} -r1 ${r1_args} -r2 ${r2_args}
    """
}
