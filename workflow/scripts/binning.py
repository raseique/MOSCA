# -*- coding: utf-8 -*-
"""
MOSCA's Binning package for clustering of 
contigs into Operational Taxonomic Units

By João Sequeira

Nov 2018
"""

from mosca_tools import run_command
import pandas as pd
import shutil
import pathlib
import snakemake


class Binner:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def run_maxbin(self, contigs, output, threads=8, reads=None, reads2=None, markerset='40', prob_threshold=0.9):
        '''
        Input:
            contigs: FASTA file with contigs
            output: basename of output
            threads: number of threads to use by Maxbin
            mg1: name of forward reads file used in the assembly
            mg2: name of reverse reads file used in the assembly
            abundance: name of abundance file (format is contig\tabundance)
            marketset: either '107' marker genes present in >95% of bacteria, or
            '40' marker gene sets that are universal among bacteria and archaea.
            '40' may be better suited for environment dominated by archaea;
            however it tends to split genomes into more bins.
        Output:
            bins named basename + .n.fasta
            abundance of each contig in each sample (mg1 and mg2) named basename + .abund1/2
            abundace of each bin for both samples named basename + .abundance
            log of workflow named basename + .log
            markergenes used to compose the bins named basename + .marker
            contigs not included in any bin named basename + .noclass
        '''
        output_folder = '/'.join(output.split('/')[:-1])  # remove the basename
        pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)
        run_command(
            f'run_MaxBin.pl -contig {contigs} -out {output} -thread {threads} -markerset {markerset} '
            f'-prob_threshold {prob_threshold}{f" -reads1 {reads}" if reads else ""}'
            f'{f" -reads2 {reads2}" if reads2 else ""}')

    def run_checkm(self, bins_folder, threads=12):
        run_command(
            f'checkm lineage_wf -x fasta -r --ali --nt -t {threads} {bins_folder} --reduced_tree '
            f'{bins_folder} --tab_table --file {bins_folder}/checkm.tsv')

    def get_bins_quality(self, table):
        table = pd.read_csv(table, sep='\t')
        hq_bins = ((table.Completeness >= 90) & (table.Contamination < 5)).sum()
        mq_bins = ((table.Completeness < 90) & (table.Completeness >= 50) & (table.Contamination < 10)).sum()
        lq_bins = ((table.Completeness < 50) & (table.Contamination < 10)).sum()
        return hq_bins, mq_bins, lq_bins

    def better_bin(self, table1, table2):
        hq_bins1, mq_bins1, lq_bins1 = self.get_bins_quality(table1)
        hq_bins2, mq_bins2, lq_bins2 = self.get_bins_quality(table2)

        if hq_bins1 > hq_bins2:
            return True
        if hq_bins1 < hq_bins2:
            return False
        if mq_bins1 > mq_bins2:
            return True
        if mq_bins1 < mq_bins2:
            return False
        if lq_bins1 > lq_bins2:
            return True
        if lq_bins1 < lq_bins2:
            return False
        return False

    def iterative_binning(self, contigs, output, threads=8, reads=None, reads2=None, markerset='40'):
        best_bin = 10
        sample = output.split('/')[-1]

        for prob_threshold in range(10, 100, 10):
            print(f'Probability threshold: {prob_threshold}')
            self.run_maxbin(contigs, f'{output}/{prob_threshold}/{sample}_{prob_threshold}', threads=threads,
                            reads=reads, reads2=reads2, markerset=markerset, prob_threshold=prob_threshold / 100)
            self.run_checkm(f'{output}/{prob_threshold}', threads=threads)

            if prob_threshold > 10:
                if self.better_bin(f'{output}/{prob_threshold}/checkm.tsv', f'{output}/{best_bin}/checkm.tsv'):
                    shutil.rmtree(f'{output}/{best_bin}')
                    print(f'Removed files for probability threshold: {best_bin} %')
                    best_bin = prob_threshold
                    print(f'New best probability threshold: {best_bin} %')
                else:
                    shutil.rmtree(f'{output}/{prob_threshold}')
                    print(f'Removed files for probability threshold: {prob_threshold} %')

        shutil.copyfile(f'{output}/{best_bin}/checkm.tsv', f'{output}/checkm.tsv')
        print(f'Best probability threshold: {best_bin} %')
        with open(f'{output}/result.txt', 'w') as f:
            f.write(f'Best probability threshold: {best_bin}')

    def run(self):
        reads = snakemake.input.reads[0]
        reads2 = snakemake.input.reads[1] if len(snakemake.input.reads) > 1 else None

        sample = snakemake.params.output.split('/')[-1]

        if snakemake.params.iterative:
            self.iterative_binning(
                snakemake.input.contigs, snakemake.params.output, threads=snakemake.threads, reads=reads, reads2=reads2,
                markerset=snakemake.params.markerset)
        else:
            self.run_maxbin(
                snakemake.input.contigs, f'{snakemake.params.output}/{sample}', threads=snakemake.threads, reads=reads,
                reads2=reads2, markerset=snakemake.params.markerset)
            self.run_checkm(snakemake.params.output, threads=snakemake.threads)


if __name__ == '__main__':
    Binner().run()
