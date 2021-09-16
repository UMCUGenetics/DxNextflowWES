# DxNextflowWES [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4551799.svg)](https://doi.org/10.5281/zenodo.4551799)

Genome Diagnostics Nextflow WES workflow

#### Get Nextflow Modules
```bash
git submodule update --init --recursive
```

#### Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

#### Running WES workflow
```bash
nextflow run WES.nf -c WES.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> --email <email> --gatk_hc_interval_list <callingtarget_dir_path> --picard_bait <baitinterval_dir_path> [-profile slurm|mac]
```

#### Running WES Fingerprint workflow
```bash
nextflow run WES_Fingerprint.nf -c WES.config --bam_path <bam_dir_path> --outdir <output_dir_path> --email <email> [-profile slurm|mac]
```

#### Create Kinship container
```bash
guixr pack -f squashfs -RR -S /bin=bin king plink-ng vcftools@0.1.14 bash glibc-utf8-locales tzdata coreutils procps grep sed bootstrap-binaries
cp change/to/guix/path.gz.squashfs kinship.gz.squashfs
unsquashfs kinship.gz.squashfs
chmod -R 0775 squashfs-root
singularity build kinship.sif squashfs-root
```
