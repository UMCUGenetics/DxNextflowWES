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