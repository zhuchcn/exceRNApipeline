import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from snakemake.shell import shell
import argparse


class HtsStats():
    def __init__(self, logfiles, names):
        self.logfiles = {
            names[i]: logfiles[i] for i in range(len(names))
        }
        self.tidy()
    
    def tidy(self):
        df = {}
        for name, logfile in self.logfiles.items():
            with open(logfile, "rt") as fh:
                raw = json.load(fh)
            data = {}
            hts_stats_keys = []
            for key in raw.keys():
                if key.startswith("hts_Stats"):
                    hts_stats_keys.append(int(key.replace("hts_Stats_", "")))
                elif key.startswith("hts_AdapterTrimmer"):
                    data["hts_AdapterTrimmer"] = raw[key]
                elif key.startswith("hts_QWindowTrim"):
                    data["hts_QWindowTrim"] = raw[key]
                elif key.startswith("hts_NTrimmer"):
                    data["hts_NTrimmer"] = raw[key]
            data["hts_Stats1"] = raw[f"hts_Stats_{min(hts_stats_keys)}"]
            data["hts_Stats2"] = raw[f"hts_Stats_{max(hts_stats_keys)}"]
            df[name] = pd.Series([
                data["hts_Stats1"]["Single_end"]["SE_in"],
                data["hts_AdapterTrimmer"]["Single_end"]["SE_out"],
                data["hts_AdapterTrimmer"]["Single_end"]["SE_discarded"],
                data["hts_AdapterTrimmer"]["Single_end"]["SE_adapterTrim"],
                data["hts_AdapterTrimmer"]["Single_end"]["SE_adapterBpTrim"],
                data["hts_QWindowTrim"]["Single_end"]["SE_out"],
                data["hts_QWindowTrim"]["Single_end"]["SE_rightTrim"],
                data["hts_QWindowTrim"]["Single_end"]["SE_leftTrim"],
                data["hts_QWindowTrim"]["Single_end"]["SE_discarded"],
                data["hts_NTrimmer"]["Single_end"]["SE_out"],
                data["hts_NTrimmer"]["Single_end"]["SE_leftTrim"],
                data["hts_NTrimmer"]["Single_end"]["SE_rightTrim"],
                data["hts_NTrimmer"]["Single_end"]["SE_discarded"],
                data["hts_Stats2"]["Single_end"]["SE_out"]
            ], index=[
                "input",
                "AdapterTrimmer_out",
                "AdapterTrimmer_discarded",
                "AdapterTrimmer_adapterTrim",
                "AdapterTrimmer_adapterBpTrim",
                "QWindowTrim_out",
                "QWindowTrim_rightTrim",
                "QWindowTrim_leftTrim",
                "QWindowTrim_discard",
                "NTrimmer_out",
                "NTrimmer_leftTrim",
                "NTrimmer_rightTrim",
                "NTrimmer_discarded",
                "output"
            ])
        self.df = pd.DataFrame(df)
    
    def barplot(self, file):
        sns.set()
        sns.set_style("ticks")
        f, ax = plt.subplots(figsize=(10, 6), dpi=200)
        self.df.loc[["output", "AdapterTrimmer_discarded", 
                          "QWindowTrim_discard", "NTrimmer_discarded"],].T\
            .plot(kind="bar", stacked=True)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(file, dpi=200, bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
    
    def to_tsv(self, path):
        self.df.to_csv(path, sep="\t")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log-files", type=str, nargs="+")
    parser.add_argument("-n", "--sample-names", type=str, nargs="+")
    parser.add_argument("-o", "--output-tsv", type=str)
    parser.add_argument("-p", "--output-plot", type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    stats = HtsStats(args.log_files, args.sample_names)
    if args.output_tsv:
        stats.to_tsv(args.output_tsv)
    if args.output_plot:
        stats.barplot(args.output_plot)


if __name__ == "__main__":
    main()