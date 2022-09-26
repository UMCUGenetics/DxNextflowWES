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