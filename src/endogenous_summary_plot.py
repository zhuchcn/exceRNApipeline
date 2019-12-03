import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class EndogenousCountSummary():
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(self.filepath, sep="\t", index_col=0)
    
    def stacked_barplot(self, filepath):
        fig = plt.figure(figsize=(10, 6), dpi = 200)
        heights = self.df.apply(lambda x: x/sum(x), axis=0)
        bottoms = heights.cumsum(axis=0)
        for index in range(heights.shape[0]):
            height = heights.iloc[index, :]
            bottom = None if index == 0 else bottoms.iloc[index-1,:]
            plt.bar(
                heights.columns,
                height=height,
                width=0.6,
                bottom=bottom,
                label=heights.index[index]
            )
        plt.xticks(heights.columns, heights.columns, rotation="vertical")
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(filepath, dpi = 200,bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')


if __name__ == "__main__":
    ecs = EndogenousCountSummary(snakemake.input[0])
    ecs.stacked_barplot(snakemake.output[0])