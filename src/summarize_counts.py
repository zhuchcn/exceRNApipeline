from itertools import zip_longest
import pandas as pd


def merge_counts(counts, output, attrs, reduce=False):
    handles = {}
    for sample, path in counts.items():
        handles[sample] = open(path, 'rt')
    with open(output, 'w') as fh:
        samples = list(counts.keys())
        headers = attrs + '\t' + '\t'.join(samples) + '\n'
        fh.write(headers)
        for lines in zip_longest(*[handles[sample] for sample in samples]):
            if lines[0].startswith('__'):
                continue
            values = [line.rstrip().split('\t')[-1] for line in lines]
            gene_annos = '\t'.join(lines[0].rstrip().split('\t')[:-1])
            if reduce:
                if len([v for v in values if int(v) > 0]) < 1:
                    continue
            line = gene_annos + '\t' + '\t'.join(values) + '\n'
            fh.write(line)
    for h in handles.values():
        h.close()

def summarize_counts(merged_counts, output, key, samples):
    summarized = None
    for path in merged_counts:
        df = pd.read_csv(
            path, sep='\t', index_col=False
        )[['gene_type'] + samples]
        # this deals with the error when the merged_count is empty.
        if df.shape[0] == 0:
            print(f"No gene count found in {path}", flush=True)
            continue
        df = df.melt(key, var_name='sample')\
            .groupby([key, 'sample'])\
            .agg({'value': sum})\
            .pivot_table(index='gene_type', columns='sample')
        df.columns = df.columns.droplevel()
        df.columns.name = None
        df.index.name = None
        if not summarized:
            summarized = df
        else:
            summarized = summarized.append(df)
    if summarized is None:
        summarized = df
        print("No gene count in any sample was saved. Please verify your data.",
              flush=True)
    summarized.to_csv(output, sep='\t')

def main():
    for count_type in snakemake.input.keys():
        counts = {}
        for sample, path in zip_longest(snakemake.params.samples,
                    list(snakemake.input.get(count_type))):
            counts[sample] = path
        output = snakemake.output.get(count_type)
        # attrs are the first several columns that are not gene counts.
        attrs = snakemake.params.get("attrs") \
                    or 'gene_id\tgene_name\tgene_type'
        merge_counts(counts, output, attrs, reduce=True)
    merged_counts = [
        snakemake.output.get(key) for key in snakemake.input.keys()
    ]
    # key is the attribute used to summarized gene counts.
    key = snakemake.params.get('key') or 'gene_type'
    samples = snakemake.params.samples
    summarize_counts(merged_counts, snakemake.output.summary, key, samples)

if __name__ == '__main__':
    main()
