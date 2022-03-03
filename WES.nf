#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

// Mapping modules
include BWAMapping from './NextflowModules/BWA-Mapping/bwa-0.7.17_samtools-1.9/Mapping.nf' params(
    genome_fasta: "$params.genome", optional: '-c 100 -M'
)
include MarkdupMerge as Sambamba_MarkdupMerge from './NextflowModules/Sambamba/0.7.0/Markdup.nf'

// IndelRealignment modules
include RealignerTargetCreator as GATK_RealignerTargetCreator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/RealignerTargetCreator.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: "$params.gatk_rtc_options"
)
include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome", optional: ""
)
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

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

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome",
    optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES"
)

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(
    genome: "$params.genome",
    optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(
    optional: "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
)
include CollectHsMetrics as PICARD_CollectHsMetrics from './NextflowModules/Picard/2.22.0/CollectHsMetrics.nf' params(
    genome: "$params.genome", bait:"$params.dxtracks_path/$params.picard_bait",
    target: "$params.dxtracks_path/$params.picard_target",
    optional: "METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include MultiQC from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(
    optional: "--config $baseDir/assets/multiqc_config.yaml"
)
include VerifyBamID2 from './NextflowModules/VerifyBamID/2.0.1--h32f71e1_2/VerifyBamID2.nf'

//SNParray calling
include IntervalListTools as PICARD_IntervalListToolsSNP from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count:"100", optional: ""
)
include HaplotypeCallerGVCF as GATK_HaplotypeCallerGVCF from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(
    genome:"$params.genome", emit_ref_confidence: "BP_RESOLUTION", compress:true, optional: ""
)
include MergeGvcfs as GATK_MergeGvcfs from './NextflowModules/GATK/4.2.1.0/MergeVcfs.nf' params(
    genome:"$params.genome", compress:true
)

//Genotype gVCF SNParray calling
include GenotypeGVCF as GATK_GenotypeGVCF from "./NextflowModules/GATK/4.2.1.0/GenotypeGvcfs.nf" params(
    genome:"$params.genome",  optional: "-all-sites", compress:true
)
include MergeVcfs as GATK_MergeVcfs from './NextflowModules/GATK/4.2.1.0/MergeVcfs.nf' params(
    genome:"$params.genome", compress:true
)

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

// Define chromosomes used to scatter GATK_RealignerTargetCreator
def chromosomes = Channel.fromPath(params.genome.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

// Define ped file, used in Kinship
def ped_file = file("${params.ped_folder}/${analysis_id}.ped")
if (!ped_file.exists()) {
    exit 1, "ERROR: ${ped_file} not found."
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

    // GATK UnifiedGenotyper (fingerprint)
    GATK_UnifiedGenotyper(Sambamba_Merge.out)

    // Clarity epp
    ClarityEppIndications(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> sample_id})

    // ExonCov
    ExonCovImportBam
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]}
    )
    ExonCovSampleQC(
        ExonCovImportBam.out.join(ClarityEppIndications.out)
            .map{sample_id, exoncov_id, indication -> [analysis_id, exoncov_id, indication]}
            .groupTuple()
    )

    // ExomeDepth
    rg_id = BWAMapping.out.map{sample_id, rg_id, bam_file, bai_file -> [sample_id, rg_id.split('_')[1]]}.unique().groupTuple()

    GetRefset(rg_id)
    ExomeDepth(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, sample_id, bam_file, bai_file]})
    ExomeDepthSummary(analysis_id, ExomeDepth.out.HC_stats_log.collect())

    // Kinship
    Kinship(GATK_CombineVariants.out)

    // QC
    FastQC(fastq_files)

    PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    PICARD_EstimateLibraryComplexity(Sambamba_Merge.out)
    PICARD_CollectHsMetrics(Sambamba_Merge.out)
    CreateHSmetricsSummary(PICARD_CollectHsMetrics.out.collect())

    Sambamba_Flagstat(Sambamba_Merge.out)
    GetStatsFromFlagstat(Sambamba_Flagstat.out.collect())

    VerifyBamID2(Sambamba_Merge.out.groupTuple())

    MultiQC(analysis_id, Channel.empty().mix(
        FastQC.out,
        PICARD_CollectMultipleMetrics.out,
        PICARD_EstimateLibraryComplexity.out,
        PICARD_CollectHsMetrics.out,
        VerifyBamID2.out.map{sample_id, self_sm -> [self_sm]},
        ExonCovSampleQC.out
    ).collect())

    // GATK HaplotypeCaller SNParray target
    PICARD_IntervalListToolsSNP(Channel.fromPath("$params.dxtracks_path/$params.gatk_hc_interval_list_snparray"))
    GATK_HaplotypeCallerGVCF(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}
        .groupTuple()
        .combine(PICARD_IntervalListToolsSNP.out.flatten())
    )
    GATK_MergeGvcfs(
        GATK_HaplotypeCallerGVCF.out.map{output_name, vcf_files, vcf_idx_files, interval_file -> [output_name, vcf_files, vcf_idx_files]}
        .groupTuple()
    )

    // Genotyping GVCF SNParray target
    GATK_GenotypeGVCF(GATK_HaplotypeCallerGVCF.out)
    GATK_MergeVcfs(GATK_GenotypeGVCF.out.groupTuple())

    // BAF analysis
    BAF_IGV(GATK_MergeVcfs.out)

    // UPD analysis
    ParseChildFromFullTrio(ped_file, Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]})
    UPD_IGV(ped_file, analysis_id, 
        ParseChildFromFullTrio.out.splitCsv().flatten(), 
        GATK_MergeVcfs.out.map{output_name, vcf_files, vcf_idx_files -> [vcf_files]}.collect()
    )

    // IGV sessions
    Single_IGV(GetRefset.out.map{sample_id, refset -> [sample_id, refset.split('\n').join(''), analysis_id]})
    Family_IGV(ParseChildFromFullTrio.out.splitCsv()
        .flatten()
        .cross(GetRefset.out.map{sample_id, refset -> [sample_id, refset.split('\n').join('')]})
        .map{sample_id, refset -> [sample_id, refset[1], analysis_id, ped_file]}
    )

    TrendAnalysisTool(
        GATK_CombineVariants.out.map{id, vcf_file, idx_file -> [id, vcf_file]}
            .concat(GetStatsFromFlagstat.out.map{file -> [analysis_id, file]})
            .concat(CreateHSmetricsSummary.out.map{file -> [analysis_id, file]})
            .groupTuple()
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

// Custom processes
process ExonCovImportBam {
    // Custom process to run ExonCov import_bam
    tag {"ExonCov ImportBam ${sample_id}"}
    label 'ExonCov'
    label 'ExonCov_ImportBam'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py import_bam \
        --threads ${task.cpus} \
        --overwrite \
        --print_sample_id \
        --exon_bed ${params.dxtracks_path}/${params.exoncov_bed} \
        ${analysis_id} WES ${bam_file} | tr -d '\n'
        """
}

process ExonCovSampleQC {
    // Custom process to run ExonCov sample_qc
    tag {"ExonCov Sample QC ${analysis_id}"}
    label 'ExonCov'
    label 'ExonCov_SampleQC'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_ids, indications)

    output:
        path("${analysis_id}.ExonCovQC_check.out")

    script:
        def samples = sample_ids.collect{"$it"}.join(" ")
        def panels = indications.collect{"$it"}.join(" ")
        """
        source ${params.exoncov_path}/venv/bin/activate
        python ${params.exoncov_path}/ExonCov.py sample_qc \
        -s ${samples} -p ${panels} > ${analysis_id}.ExonCovQC_check.out
        """
}

process ClarityEppIndications {
    // Custom process to run clarity_epp export sample_indications
    tag {"ClarityEppExportSampleIndications ${analysis_id}"}
    label 'ClarityEpp'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a clarity export restarting the workflow.

    input:
        val(sample_id)

    output:
        tuple(sample_id, stdout)

    script:
        """
        source ${params.clarity_epp_path}/venv/bin/activate
        python ${params.clarity_epp_path}/clarity_epp.py export sample_indications \
        -a ${sample_id} | cut -f 2 | grep -v 'Indication' | tr -d '\n'
        """
}

process GetRefset{
    // Custom process to run get exomedepth reference set from exomedepth db
    tag {"GetRefset ${sample_id}"}
    label 'GetRefset'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false

    input:
        tuple(sample_id, rg_id)

    output:
        tuple(sample_id, stdout)

    script:
        def rg_ids = rg_id.collect().join(" ")
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/exomedepth_db.py add_sample_return_refset ${sample_id} ${rg_ids}
        """
}

process ExomeDepth {
    // Custom process to run Exomedepth
    tag {"ExomeDepth ${sample_id}"}
    label 'ExomeDepth'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, sample_id, path(bam_file), path(bai_file))

    output:
        path("*.log", emit: ED_log)
        path("HC_*.igv", emit: HC_igv)
        path("UMCU_*.igv", emit: UMCU_igv)
        path("HC_*.vcf", emit: HC_vcf)
        path("UMCU_*.vcf", emit: UMCU_vcf)
        path('HC_*_stats.log', emit: HC_stats_log)
        path('UMCU_*_stats.log', emit: UMCU_stats_log)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/run_ExomeDepth.py callcnv ./ ${bam_file} ${analysis_id} ${sample_id}
        """
}

process ExomeDepthSummary {
    // Custom process to stats from ExomeDepth analysis
    tag {"ExomeDepthSummary"}
    label 'ExomeDepthSummary'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(analysis_id)
        path(exomedepth_logs)

    output:
        path("${analysis_id}_exomedepth_summary.txt")

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/exomedepth_summary.py ${exomedepth_logs}  > ${analysis_id}_exomedepth_summary.txt
        """
}


process Kinship {
    // Custom process to run Kinship tools
    // Container does not work
    // king: run.c:355: main: Unexpected error: No such file or directory.
    // Aborted (core dumped)
    tag {"Kinship ${analysis_id}"}
    label 'Kinship'
    //container = '/hpc/diaggen/software/guix_containers/kinship.sif'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(analysis_id, path(vcf_file), path(vcf_index))

    output:
        tuple(analysis_id, path("${analysis_id}.kinship"), path("${analysis_id}.kinship_check.out"))

    script:
        """
        ${params.vcftools_path}/vcftools --vcf ${vcf_file} --plink
        ${params.plink_path}/plink --file out --make-bed --noweb
        ${params.king_path}/king -b plink.bed --kinship
        cp king.kin0 ${analysis_id}.kinship
        python ${baseDir}/assets/check_kinship.py ${analysis_id}.kinship ${ped_file} > ${analysis_id}.kinship_check.out
        """
}

process GetStatsFromFlagstat {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"GetStatsFromFlagstat"}
    label 'GetStatsFromFlagstat'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(flagstat_files)

    output:
        path('run_stats.txt')

    script:
        """
        python ${baseDir}/assets/get_stats_from_flagstat.py ${flagstat_files} > run_stats.txt
        """
}

process CreateHSmetricsSummary {
    // Custom process to run get_stats_from_flagstat.pl
    tag {"CreateHSmetricsSummary"}
    label 'CreateHSmetricsSummary'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hsmetrics_files)

    output:
        path('HSMetrics_summary.txt')

    script:
        """
        python ${baseDir}/assets/create_hsmetrics_summary.py ${hsmetrics_files} > HSMetrics_summary.txt
        """
}
process BAF_IGV {
    // Custom process to run BAF analysis
    tag {"BAF_IGV ${output_name}"}
    label 'BAF_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(output_name, path(vcf_files), path(vcf_idx_files))

    output:
        path("${output_name}_baf.igv")

    script:
        """
        source ${params.baf_path}/venv/bin/activate
        python ${params.baf_path}/make_BAF_igv.py ${vcf_files} -o ${output_name}_baf.igv -c
        """
}

process ParseChildFromFullTrio {
    //Custom process to parse PED file and output sampleID of children with both parents.
    tag {"ParseChildFromFullTrio ${analysis_id}"}
    label 'ParseChildFromFullTrio'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(ped_file)
        val(vcf_files)

    output:
        stdout emit: trio_sample

    script:
        """
        python ${baseDir}/assets/parse_child_from_fulltrio.py ${ped_file} ${vcf_files}
        """
}

process UPD_IGV {
    // Custom process to run UPD analysis
    tag {"UPD_IGV $trio_sample"}
    label 'UPD_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(ped_file)
        val(analysis_id)
        val(trio_sample)
        path(vcf_files)

    output:
        path("*.igv", emit: UPD_IGV_files)

    script:
        """
        source ${params.upd_path}/venv/bin/activate
        python ${params.upd_path}/make_UPD_igv.py ${ped_file} ${analysis_id} $trio_sample ${vcf_files} -c
        """
}

process Single_IGV {
    // Custom process to run Single sample IGV analysis
    tag {"Single_IGV $sample_id"}
    label 'Single_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, refset, analysis_id)

    output:
        path("*.xml", emit: Single_IGV_file)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/igv_xml_session.py single_igv ./ ${sample_id} ${analysis_id} ${refset}
        """
}


process Family_IGV {
    // Custom process to run Family IGV analysis
    tag {"Family_IGV $sample_id"}
    label 'Family_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, refset, analysis_id, ped_file)

    output:
        path("*.xml", emit: Family_IGV_file)

    script:
        """
        source ${params.exomedepth_path}/venv/bin/activate
        python ${params.exomedepth_path}/igv_xml_session.py family_igv ./ ${ped_file} ${analysis_id} ${sample_id} ${refset}
        """
}

process TrendAnalysisTool {
    // Custom process to run Trend_Analysis_tool
    tag {"TrendAnalysisTool ${analysis_id}"}
    label 'TrendAnalysisTool'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(analysis_id, path(input_files))

    script:
        """
        source ${params.trend_analysis_path}/venv/bin/activate
        python ${params.trend_analysis_path}/trend_analysis.py upload processed_data ${analysis_id} .
        """
}

process SavePedFile {
    tag {"SavePedFile ${analysis_id}"}
    label 'SavePedFile'
    shell = ['/bin/bash', '-euo', 'pipefail']
    cache = false  //Disable cache to force a new ped file copy when restarting the workflow.

    output:
        path("*.ped")

    script:
        """
        cp ${ped_file} ./
        """
}

process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowWes' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_tracks' >> repository_version.log
        git --git-dir=${params.dxtracks_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'ExonCov' >> repository_version.log
        git --git-dir=${params.exoncov_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'clarity_epp' >> repository_version.log
        git --git-dir=${params.clarity_epp_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'ExomeDepth' >> repository_version.log
        git --git-dir=${params.exomedepth_path}/../.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_UPD' >> repository_version.log
        git --git-dir=${params.upd_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_BAF' >> repository_version.log
        git --git-dir=${params.baf_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'TrendAnalysis' >> repository_version.log
        git --git-dir=${params.trend_analysis_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
