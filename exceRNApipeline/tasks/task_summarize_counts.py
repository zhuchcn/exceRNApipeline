from itertools import zip_longest
import pandas as pd
import argparse


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
        if summarized is None:
            summarized = df
        else:
            summarized = summarized.append(df)
    if summarized is None:
        summarized = df
        print("No gene count in any sample was saved. Please verify your data.",
              flush=True)
    summarized.to_csv(output, sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-gencode', type=str, nargs="+")
    parser.add_argument('--input-tRNA', type=str, nargs="+")
    parser.add_argument('--input-piRNA', type=str, nargs="+")
    parser.add_argument('--output-gencode', type=str)
    parser.add_argument('--output-tRNA', type=str)
    parser.add_argument('--output-piRNA', type=str)
    parser.add_argument('--output-summary', type=str)
    parser.add_argument('--sample-names', type=str, nargs="+")
    parser.add_argument('--output-attrs', type=str, nargs='*')
    parser.add_argument('--key', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    input_keys = [key.replace('input_', '') for key in args.__dict__.keys() 
                    if key.startswith('input_')]
    for count_type in input_keys:
        counts = {}
        for sample, path in zip_longest(args.sample_names,
                    list(args.__dict__[f'input_{count_type}'])):
            counts[sample] = path
        output = args.__dict__[f'output_{count_type}']
        # attrs are the first several columns that are not gene counts.
        attrs = args.output_attrs \
                    or 'gene_id\tgene_name\tgene_type'
        merge_counts(counts, output, attrs, reduce=True)
    merged_counts = [args.__dict__[f'output_{key}'] for key in input_keys]
    # key is the attribute used to summarized gene counts.
    key = args.key or 'gene_type'
    samples = args.sample_names
    summarize_counts(merged_counts, args.output_summary, key, samples)

if __name__ == '__main__':
    main()
