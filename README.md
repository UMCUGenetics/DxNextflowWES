# DxNextflowWorkflows
Genome Diagnostics Nextflow workflows

#### Get Nextflow Modules
```bash
git submodule update --init --recursive
```

#### Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

#### Running fingerprint pipeline
```bash
tools/nextflow run MipFingerprint.nf -c MipFingerprint.config --samplesheet samples.tsv --outdir <output_dir_name> [-profile slurm|sge]
```

#### Running WES workflow
```bash
tools/nextflow run WES.nf -c WES.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> [-profile slurm|sge]
```
