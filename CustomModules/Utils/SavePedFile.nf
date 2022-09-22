process SavePedFile {
    tag {"SavePedFile ${analysis_id}"}
    label 'SavePedFile'
    shell = ['/bin/bash', '-euo', 'pipefail']
    cache = false  //Disable cache to force a new ped file copy when restarting the workflow.

    input:
        path(ped_file)

    output:
        path("*.ped", includeInputs=True)

    script:
        """
        cp --remove-destination "\$(readlink ${ped_file})" ./
        """
}