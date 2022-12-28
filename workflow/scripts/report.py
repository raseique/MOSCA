# -*- coding: utf-8 -*-
"""
General report construction and export

By João Sequeira

Oct 2019
"""

import argparse
import glob
import pandas as pd
from mosca_tools import parse_blast, run_pipe_command, parse_fastqc_report, count_on_file
from zipfile import ZipFile


class Reporter:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's technical and quality control reports")
        parser.add_argument("-e", "--experiments", required=True, help="Experiments file")
        parser.add_argument("-o", "--output", help="Output directory")
        parser.add_argument("-ldir", "--lists-directory", help="Directory with lists for Reporter")
        parser.add_argument("-if", "--input-format", default='tsv', choices=['tsv', 'excel'])
        parser.add_argument("--no-differential-expression", action='store_true', default=False)
        parser.add_argument("-s", "--suffix", default="",
                            help="If files don't end with _R1, _R2, this will be added to ends")
        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    def write_technical_report(self, output):
        """
        Input:
            output: str - filename to write tools and respective versions
        Output:
            a file named [output] will be written with information concerning
            the softwares used by MOSCA and respective versions
        """
        conda_list = run_pipe_command('conda list', output='PIPE').split('\n')[2:]
        lines = [line.split() for line in conda_list]
        lines[0] = lines[0][1:]
        df = pd.DataFrame(lines, columns=lines.pop(0)).set_index('Name')
        df[['Version']].to_csv(output, sep='\t')

    def initialize_report(self, reporter_columns):
        print('Initializing Report')
        self.report = pd.DataFrame(columns=reporter_columns)

    def info_from_fastqc(self, output_dir, name, col_name, prefix, prefix2terms, fastq_columns):
        """
        Input:
            name: str - name of MG/MT sample
            prefix: str - [Initial quality assessment], [Before quality trimming] or
            [After quality trimming] - the step that concerns the FastQC report
            performed_rrna_removal: bool - True if rRNA removal was performed, False
            otherwise. Should be False for MG, and True for MT
        Output:
            self.report will be updated with information from the report, on the line
            named 'name', and the columns that start with 'prefix'
        """
        reports = [parse_fastqc_report(f'{output_dir}/Preprocess/FastQC/{prefix2terms[prefix][0]}{name}_'
                                       f'{prefix2terms[prefix][i]}_fastqc/fastqc_data.txt') for i in [1, 2]]
        self.report.loc[col_name, f'{prefix} # of reads'] = reports[0]['Basic Statistics'][1].loc[
            'Total Sequences']['Value']

        for column in fastq_columns:
            for i in range(2):
                if column not in reports[i].keys():
                    reports[i][column] = (
                        'Not available', None)  # only the not available matters. And nothing else matters!...
            if reports[0][column][0] == reports[1][column][0]:
                self.report.loc[col_name, f'{prefix} {column}'] = reports[0][column][0]
            else:
                self.report.loc[col_name, f'{prefix} {column}'] = (
                    f'{reports[0][column][0]} (forward) {reports[1][column][0]} (reverse)')

    def info_from_preprocessing(self, output_dir, name, input_file, fastq_columns, performed_rrna_removal=False,
                                suffix=''):
        print(f'Retrieving preprocessing information for dataset: {name}')
        if name not in self.report.index:
            # add a new line to the report, with the name of the dataset
            self.report.loc[name] = [None] * len(self.report.columns)
        self.report.loc[name] = self.report.loc[name].fillna(value='')

        adapter_files = open(f'{output_dir}/Preprocess/Trimmomatic/{name}_adapters.txt').read().split('\n')
        if len(adapter_files[0]) > 0 and not adapter_files[0] in ['None', 'Fail']:
            adapter = adapter_files[0].split('/')[-1].split('.fa')[0]
        else:
            adapter_files = []
            adapter = None

        # For each preprocessing step, a tuple of (prefix, suffix for forward, suffix for reverse)
        prefix2terms = {'[Initial quality assessment]': ('', f'R1{suffix}', f'R2{suffix}'),
                        '[Before quality trimming]': (
                            ('norrna_', 'fwd', 'rev') if performed_rrna_removal else (
                            'noadapters_', f'{adapter}_forward_paired', f'{adapter}_reverse_paired')
                        if adapter is not None else ('', 'R1', 'R2')),
                        '[After quality trimming]': ('quality_trimmed_', 'forward_paired', 'reverse_paired')}

        # Initial assessment
        if '_R' in input_file:  # is paired-end
            ends_name = input_file.split('/')[-1].split('_R')[0]
        else:  # is single-end
            ends_name = input_file.split('/')[-1].split('.f')[0]

        self.info_from_fastqc(output_dir, ends_name, name, '[Initial quality assessment]', prefix2terms, fastq_columns)

        # After adapter removal
        try:
            if len(adapter_files) > 0:
                self.report.loc[name, '[Adapter removal] adapter files'] = ', '.join(set(adapter_files))
            else:
                self.report.loc[name, '[Adapter removal] adapter files'] = 'None'
        except:
            print('Failed at adapter removal!')
            self.report.to_csv(f'{output_dir}/report.tsv', sep='\t')

        # Quality trimming
        if adapter is not None or performed_rrna_removal:
            right_name = name
        else:  # still using original files, no rRNA or adapter removal happened
            right_name = ends_name

        self.info_from_fastqc(output_dir, right_name, name, '[Before quality trimming]', prefix2terms, fastq_columns)

        self.report.loc[name, '[Quality trimming] Parameters'] = '; '.join([
            file for file in set(open(f'{output_dir}/Preprocess/Trimmomatic/{name}_quality_params.txt').read(
            ).split('\n')) if
            len(file) > 2])  # TODO - because '' must be interfering, try to cut the problem at the root before troubles

        self.info_from_fastqc(output_dir, name, name, '[After quality trimming]', prefix2terms, fastq_columns)

    def set_samples(self, experiments):
        self.report = pd.merge(experiments[['Name', 'Sample']], self.report, left_on='Name', right_index=True,
                               how='outer')

    def info_from_assembly(self, output_dir, sample):
        print(f'Retrieving assembly information for sample: {sample}')
        qc_report = pd.read_csv(f'{output_dir}/Assembly/{sample}/quality_control/report.tsv', sep='\t', index_col=0
                                ).transpose()
        qc_report.index = [sample]

        for col in qc_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, f'[Assembly] {col}'] = qc_report.loc[sample][col]

    def info_from_annotation(self, output_dir, sample):
        print(f'Retrieving annotation information for sample: {sample}')
        sample_report = dict()
        sample_report['# of proteins detected'] = count_on_file('>', f'{output_dir}/Annotation/{sample}/fgs.faa')

        sample_report['# of proteins annotated (DIAMOND)'] = (
                len(set(parse_blast(f'{output_dir}/Annotation/{sample}/aligned.blast')['qseqid'])) - (
            count_on_file('*', f'{output_dir}/Annotation/{sample}/aligned.blast')))

        sample_report['# of proteins annotated (reCOGnizer)'] = (
            len(set(pd.read_csv(f'{output_dir}/Annotation/{sample}/COG_report.tsv', sep='\t')['qseqid'])))

        sample_report = pd.DataFrame.from_dict(sample_report, orient='index').transpose()
        sample_report.index = [sample]

        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, f'[Annotation] {col}'] = (
                sample_report.loc[sample][col])

    def info_from_binning(self, output_dir, sample):
        print(f'Retrieving binning information for sample: {sample}')
        sample_report = dict()
        sample_report['# of bins'] = len(glob.glob(f'{output_dir}/Binning/{sample}/{sample}.*.fasta'))
        checkm = pd.read_csv(f'{output_dir}/Binning/{sample}/checkm.tsv', sep='\t')
        sample_report['# of high-quality drafts'] = (
                (checkm['Completeness'] > 90) & (checkm['Contamination'] < 5)).sum()
        sample_report['# of medium-quality drafts'] = (
                (checkm['Completeness'] < 90) & (checkm['Completeness'] > 50) & (checkm['Contamination'] < 10)).sum()
        sample_report['# of low-quality drafts'] = (
                (checkm['Completeness'] < 50) & (checkm['Contamination'] < 10)).sum()
        sample_report = pd.DataFrame.from_dict(sample_report, orient='index').transpose()
        sample_report.index = [sample]

        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, f'[Binning] {col}'] = (
                sample_report.loc[sample][col])

    def info_from_alignment(self, p_report, mt_name):
        self.report.set_index('Name', inplace=True)
        self.report.loc[mt_name, '[Metatranscriptomics] # of reads aligned'] = p_report[mt_name].sum()
        self.report.reset_index(inplace=True)

    def info_from_differential_expression(self, output_dir, sample, cutoff=0.01):
        de_results = pd.read_csv(
            f'{output_dir}/Quantification/{sample}/condition_treated_results.tsv', index_col=0, sep='\t')
        self.report.loc[self.report['Sample'] == sample, '[Gene expression] # of differentially expressed proteins'] = (
            (de_results['padj'] < cutoff).sum())

    def zip_files(self, files, output):
        with ZipFile(output, 'w') as archive:
            for file in files:
                archive.write(file)

    def run(self):
        args = self.get_arguments()

        fastq_columns = open(f'{args.lists_directory}/fastqc_columns.txt').read().split('\n')
        reporter_columns = open(f'{args.lists_directory}/reporter_columns.txt').read().split('\n')

        self.write_technical_report(f'{args.output}/technical_report.tsv')

        exps = (
            pd.read_csv(args.experiments, sep='\t') if args.input_format == 'tsv' else pd.read_excel(args.experiments))

        self.initialize_report(reporter_columns)

        for i in exps.index:
            self.info_from_preprocessing(
                args.output, exps.iloc[i]['Name'], exps.iloc[i]['Files'].split(',')[0], fastq_columns,
                performed_rrna_removal=(False if exps.iloc[i]['Data type'] == 'dna' else True), suffix=args.suffix)

        self.set_samples(exps)

        for sample in set(exps['Sample']):

            self.info_from_assembly(args.output, sample)

            self.info_from_annotation(args.output, sample)

            self.info_from_binning(args.output, sample)

            if not args.no_differential_expression:
                self.info_from_differential_expression(args.output, sample)

            with open(f"{args.output}/{sample}_report.txt", 'w') as f:
                f.write('done')

        protein_report = pd.read_excel(f'{args.output}/MOSCA_Protein_Report.xlsx')
        for mt_name in exps[exps["Data type"] == 'mrna']['Name']:
            self.info_from_alignment(protein_report, mt_name)

        self.report.to_excel(f'{args.output}/MOSCA_General_Report.xlsx', index=False)

        self.zip_files([f'{args.output}/{filename}' for filename in [
                'MOSCA_Protein_Report.xlsx',
                'MOSCA_Entry_Report.xlsx',
                'MOSCA_General_Report.xlsx',
                'technical_report.tsv'
            ]],
            f'{args.output}/MOSCA_results.zip')


if __name__ == '__main__':
    Reporter().run()
