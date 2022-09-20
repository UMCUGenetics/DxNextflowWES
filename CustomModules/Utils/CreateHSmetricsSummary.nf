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