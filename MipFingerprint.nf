#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include 'NextflowModules/Utils/fastq.nf'
include FastQC from 'NextflowModules/FastQC/0.11.8/FastQC.nf' params(params)
include MipsTrimDedup from 'NextflowModules/Mips/1.0.1/MipsTrimDedup.nf' params(
    outdir: params.outdir,
    mips_trim_dedup: params.mips_trim_dedup,
    design_file: params.mip_design_file,
    uuid_length: params.mip_uuid_length,
    uuid_read: params.mip_uuid_read,
    )
include MEM as BWAMEM from 'NextflowModules/BWA/0.7.17/MEM.nf' params(params)
include view_sort as Sambamba_view_sort from 'NextflowModules/Sambamba/0.7.0/view_sort.nf' params(params)

fastq_files = extractFastqFromDir(params.fastq_path)
samples = fastq_files.groupTuple(by:[0])

FastQC(fastq_files)
MipsTrimDedup(samples)
BWAMEM(MipsTrimDedup.out)
SambambaViewSort(BWAMEM.out)
