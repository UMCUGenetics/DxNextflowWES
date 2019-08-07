#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include MipsTrimDedup from 'NextflowModules/Mips/1.0.1/MipsTrimDedup.nf' params(
    outdir: params.outdir,
    mips_trim_dedup: params.mips_trim_dedup,
    design_file: params.mip_design_file,
    uuid_length: params.mip_uuid_length,
    uuid_read: params.mip_uuid_read,
    )
include FastQC from 'NextflowModules/FastQC/0.11.8/FastQC.nf' params(params)
include MEM as BWA_MEM from 'NextflowModules/BWA/0.7.17/MEM.nf' params(params)
include UnifiedGenotyper as GATK_UnifiedGenotyper_fingerprint from 'NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(
    outdir: params.outdir,
    process_outdir: 'fingerprint',
    gatk: params.gatk,
    genome: params.genome,
    intervals: params.fingerprint_target,
    dbsnp: params.fingerprint_target,
    output_mode: 'EMIT_ALL_SITES'
    )

params.samplesheet
params.outdir

Channel.fromPath( file(params.samplesheet) )
    .splitCsv(header: true, sep: '\t')
    .map{row ->
        def sample = row['Sample']
        def reads1 = row['R1'].tokenize( ',' ).collect { file(it) }
        def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
        return [ sample, reads1, reads2 ]
    }
    .tap{samples_R1_R2_fastq} // set of all fastq R1 R2 per sample
    .map { sample_ID, reads1, reads2 ->
        return [ sample_ID, [reads1, reads2 ].flatten() ]
    }
    .tap{samples_all_fastq} // set of all fastq per sample

MipsTrimDedup(samples_R1_R2_fastq)
FastQC(samples_all_fastq)
BWA_MEM(mips_trim_dedup.out)
GATK_UnifiedGenotyper_fingerprint(bwa_mem.out)
