# DxNextflowWES
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
tools/nextflow run WES.nf -c WES.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> [-profile slurm|mac]
```

#### Create Kinship container
```bash
guixr pack -f squashfs -RR -S /bin=bin king plink-ng vcftools@0.1.14 bash glibc-utf8-locales tzdata coreutils procps grep sed bootstrap-binaries
cp change/to/guix/path.gz.squashfs kinship.gz.squashfs
unsquashfs kinship.gz.squashfs
chmod -R 0775 squashfs-root
singularity build kinship.sif squashfs-root
```
