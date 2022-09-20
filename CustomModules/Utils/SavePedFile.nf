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