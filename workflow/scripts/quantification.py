# -*- coding: utf-8 -*-
"""
MOSCA's script for quantifying reads

By João Sequeira

Dec 2022
"""

import pandas as pd
import argparse
import multiprocessing
from pathlib import Path
from mosca_tools import perform_alignment, normalize_counts_by_size


def get_arguments():
    parser = argparse.ArgumentParser(description="MOSCA quantification")

    parser.add_argument("-o", "--output", help="Output directory (and input!).")
    parser.add_argument("-e", "--experiments", help="Filename of exps.")
    parser.add_argument(
        "-t", "--threads", default=multiprocessing.cpu_count() - 2, help="Number of threads to use [max available - 2]")
    args = parser.parse_args()
    args.output = args.output.rstrip('/')
    return args


def run():
    args = get_arguments()

    exps = pd.read_csv(args.experiments, sep='\t')

    for sample in set(exps['Sample']):
        mg_result = pd.DataFrame(columns=['Contig'])
        mt_result = pd.DataFrame(columns=['Gene'])
        pexps = exps[(exps['Sample'] == sample)]
        for i in pexps.index:
            if pexps.iloc[i]['Data type'] == 'mrna':
                reference = f"{args.output}/Annotation/{pexps.iloc[i]['Sample']}/fgs.ffn"
            elif pexps.iloc[i]['Data type'] == 'dna':
                reference = f"{args.output}/Assembly/{pexps.iloc[i]['Sample']}/contigs.fasta"
            else:
                continue
            if ',' in pexps.iloc[i]['Files']:
                reads = [f"{args.output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.iloc[i]['Name']}_{fr}_paired.fq"
                         for fr in ['forward', 'reverse']]
            else:
                reads = [f"{args.output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.iloc[i]['Name']}.fq"]
            perform_alignment(
                reference, reads, f"{args.output}/Quantification/{pexps.iloc[i]['Name']}", threads=args.threads)
            normalize_counts_by_size(f"{args.output}/Quantification/{pexps.iloc[i]['Name']}.readcounts", reference)
            # Read the results of alignment and add them to the readcounts result at sample level
            counts = pd.read_csv(
                f"{args.output}/Quantification/{pexps.iloc[i]['Name']}_normalized.readcounts", sep='\t',
                names=['Gene' if pexps.iloc[i]['Data type'] == 'mrna' else 'Contig', pexps.iloc[i]['Name']])
            if pexps.iloc[i]['Data type'] == 'dna':
                mg_result = pd.merge(mg_result, counts, how='outer', on='Contig')
            else:
                mt_result = pd.merge(mt_result, counts, how='outer', on='Gene')
        Path(f"{args.output}/Quantification/{sample}").mkdir(parents=True, exist_ok=True)
        if len(mg_result) > 0:
            mg_result.to_csv(f"{args.output}/Quantification/{sample}/mg.readcounts", sep='\t', index=False)
        if len(mt_result) > 0:
            mt_result.astype(int, errors='ignore').to_csv(
                f"{args.output}/Quantification/{sample}/mt.readcounts", sep='\t', index=False)


if __name__ == '__main__':
    run()
