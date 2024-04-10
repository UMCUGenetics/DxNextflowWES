#!/bin/bash

# Set parameters
java_download_url='https://download.java.net/java/GA/jdk21.0.2/f2283984656d49d69e91c558476027ac/13/GPL/openjdk-21.0.2_linux-x64_bin.tar.gz'
nextflow_version=23.10.0

# Remember the cwd
repo_dir=$PWD

echo "Updating Git submodules"

# Update submodules
git submodule update --init --recursive

# Create tools directories
mkdir -p ${repo_dir}/tools/java
mkdir -p ${repo_dir}/tools/nextflow

cd ${repo_dir}/tools/java

echo "Downloading and installing OpenJDK"
# Download OpenJDK and untar
wget ${java_download_url} -O - | tar -xz

# Only one file should be there, but still a for loop because the last part of the file is unknown
for file in ${repo_dir}/tools/java/*; do
    # Use sed to remove trailing bash, "-" and numbers
    new_name=${file%-*}
    echo "renaming" ${file} "to" ${new_name}
    # Rename the file
    mv -- ${file} ${new_name}
done

# Log the downloaded OpenJDK version + arch
java_filename=$(basename ${java_download_url})
java_base=${java_filename%%.*}
echo ${java_base} > java_version.txt

chmod +x ${repo_dir}/tools/java/jdk

${repo_dir}/tools/java/jdk/bin/java -version

echo "OpenJDK ${java_base} installed to ${repo_dir}/tools/java/jdk"

echo "Downloading and installing Nextflow"

cd ${repo_dir}/tools/nextflow

export NXF_JAVA_HOME=${repo_dir}'/tools/java/jdk'
export NXF_VER=${nextflow_version}

# Download and install nextflow
curl -s https://get.nextflow.io | bash

chmod +x ${repo_dir}/tools/nextflow/nextflow

${repo_dir}/tools/nextflow/nextflow -v

echo "Nextflow ${nextflow_version} installed on ${repo_dir}/tools/nextflow/nextflow"

echo "All done!"