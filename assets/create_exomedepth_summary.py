#! /usr/bin/env python
import argparse
import statistics

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ExomeDepth Metric summary file')
    parser.add_argument('exomedepth_logs', type=argparse.FileType('r'), nargs='*', help='Exomedepth log files')
    arguments = parser.parse_args()
    stats_dic = {"CR":[],"DP":[],"TC":[]}
    for exomedepth_qc_file in arguments.exomedepth_logs:
        for line in exomedepth_qc_file: 
            splitline = line.split()
            correlation = float(splitline[2])
            deldupratio = float(splitline[3])
            totalcount = float(splitline[4])

            stats_dic["CR"]+=[correlation]
            stats_dic["DP"]+=[deldupratio]
            stats_dic["TC"]+=[totalcount]
            print("{sample};MODEL={model};CR={correl};DP={deldupratio};TC={totalcount}\r".format(
                 sample=splitline[0],
                 model=splitline[1],
                 correl="%.4f" % correlation,
                 deldupratio="%.2f" % deldupratio,
                 totalcount="%.2f" % totalcount
                 )),
     
    print("\r")
    print("#Average_CR={}\r".format("%.4f" % statistics.mean(stats_dic["CR"]))),
    print("#Average_DP={}\r".format("%.2f" % statistics.mean(stats_dic["DP"]))),
    print("#Average_TC={}\r".format("%.2f" % statistics.mean(stats_dic["TC"]))),
    print("\r")
    print("#Median_CR={}\r".format("%.4f" % statistics.median(stats_dic["CR"]))),
    print("#Median_DP={}\r".format("%.2f" % statistics.median(stats_dic["DP"]))),
    print("#Median_TC={}\r".format("%.2f" % statistics.median(stats_dic["TC"]))),
