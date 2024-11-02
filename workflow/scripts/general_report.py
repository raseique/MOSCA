"""
MOSCA's script for producing General report

By João Sequeira

Dec 2022
"""

from mosca_tools import timed_message, blast_cols
import pandas as pd
import numpy as np
import dask.dataframe as dd
import shutil
from concurrent.futures import ThreadPoolExecutor

functional_columns = [
    'General functional category', 'Functional category', 'Protein description', 'COG ID', 'EC number (reCOGnizer)']


def make_general_report(out, exps, sample, mg_preport, mt_preport, mp_preport, de_input, did_assembly=True):
    timed_message(f'Joining data for sample: {sample}.')
    print('Reading gene calling headers.')
    with open(f'{out}/Annotation/{sample}/fgs.faa') as f:
        headers = [line.strip()[1:] for line in f if line.startswith(">")]

    report = pd.DataFrame(headers, columns=["qseqid"])

    print('Reading reCOGnizer results.')
    cog_report = dd.read_csv(f'{out}/Annotation/{sample}/COG_report.tsv', sep='\t', dtype=str)
    cog_report = cog_report[cog_report['DB ID'].str.startswith('COG') == True].rename(columns={'DB ID': 'COG ID'})
    cog_report = cog_report.groupby('qseqid').first().compute()
    report = pd.merge(report, cog_report, left_on='qseqid', right_index=True, how='left')

    print('Reading UPIMAPI results.')
    upimapi_results = dd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t', dtype=str).compute()
    report = pd.merge(upimapi_results, report, on='qseqid', how='outer')

    print('Formatting names of columns.')
    rename_cols = blast_cols + ['EC number']
    report = report.rename(columns={**{f'{col}_x': f'{col} (UPIMAPI)' for col in rename_cols},
                                    **{f'{col}_y': f'{col} (reCOGnizer)' for col in rename_cols}})

    if did_assembly:
        report['Contig'] = report['qseqid'].apply(lambda x: x.split('_')[1])

    mg_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'dna')]['Name'].tolist()
    mt_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'mrna')]['Name'].tolist()
    mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()

    def process_quantification(type_name, did_assembly, names, out, sample):
        """
        Read counts depending on the type of data (mg, mt, mp).
        """
        counts = None
        if type_name == 'mg':
            filepath = f'{out}/Quantification/{sample}_mg_norm.tsv' if did_assembly else f'{out}/Quantification/{sample}_mg.readcounts'
            counts = pd.read_csv(filepath, sep='\t', names=['Contig' if did_assembly else 'qseqid'] + names, skiprows=1)
            if did_assembly:
                counts['Contig'] = counts['Contig'].apply(lambda x: x.split('_')[1])
        elif type_name == 'mt':
            norm_filepath = f'{out}/Quantification/{sample}_mt_norm.tsv' if did_assembly else f'{out}/Quantification/{sample}_mt.readcounts'
            counts = pd.read_csv(norm_filepath, sep='\t', names=['qseqid'] + names)
        elif type_name == 'mp':
            counts = pd.read_csv(f'{out}/Metaproteomics/{sample}_mp.spectracounts', sep='\t')
            counts.rename(columns={'Main Accession': 'qseqid'}, inplace=True)
        return counts

    with ThreadPoolExecutor() as executor:
        futures = []
        if mg_names:
            futures.append(executor.submit(process_quantification, 'mg', did_assembly, mg_names, out, sample))
        if mt_names:
            futures.append(executor.submit(process_quantification, 'mt', did_assembly, mt_names, out, sample))
        if mp_names:
            futures.append(executor.submit(process_quantification, 'mp', did_assembly, mp_names, out, sample))
        results = [future.result() for future in futures]

    for result in results:
        if result is not None:
            if 'Contig' in result.columns:
                report = pd.merge(report, result, on='Contig', how='left')
            else:
                report = pd.merge(report, result, on='qseqid', how='left')

    report[mg_names + mt_names + mp_names] = report[mg_names + mt_names + mp_names].fillna(0).astype(float).astype(int)     # astype(float).astype(int) avoids "ValueError: invalid literal for int() with base 10: '2.0'"
    report.to_csv(f'{out}/MOSCA_{sample}_General_Report.tsv', sep='\t', index=False)
    return report, mg_preport, mt_preport, mp_preport, de_input


def make_general_reports(out, exps, max_lines=1000000, did_assembly=True):
    mg_report = mt_report = mp_report = de_input = pd.DataFrame(columns=['Entry'])
    writer = pd.ExcelWriter(f'{out}/MOSCA_General_Report.xlsx', engine='xlsxwriter')

    for sample in set(exps['Sample']):
        report, mg_report, mt_report, mp_report, de_input = make_general_report(
            out, exps, sample, mg_report, mt_report, mp_report, de_input, did_assembly=did_assembly)
        timed_message(f'Writing General Report for sample: {sample}.')
        if len(report) < max_lines:
            report.to_excel(writer, sheet_name=sample, index=False)
        else:
            for k, chunk in enumerate(np.array_split(report, len(report) // max_lines)):
                chunk.to_excel(writer, sheet_name=f'{sample} ({k + 1})', index=False)
    writer.close()

    timed_message('Writing quantification matrices.')
    if not mg_report.empty:
        mg_report.iloc[:, 1:] = mg_report.iloc[:, 1:].astype(float)
        mg_report = mg_report.groupby('Entry').sum().reset_index()
        mg_report.to_csv(f'{out}/Quantification/mg_entry_quant.tsv', sep='\t', index=False)
    if not mt_report.empty:
        mt_report.iloc[:, 1:] = mt_report.iloc[:, 1:].astype(float)
        mt_report = mt_report.groupby('Entry').sum().reset_index()
        mt_report.to_csv(f'{out}/Quantification/mt_entry_quant.tsv', sep='\t', index=False)
        de_input.to_csv(f'{out}/Quantification/dea_input.tsv', sep='\t', index=False)
    if not mp_report.empty:
        mp_report.iloc[:, 1:] = mp_report.iloc[:, 1:].astype(float)
        mp_report = mp_report.groupby('Entry').sum().reset_index()
        mp_report = mp_report.drop_duplicates().dropna(subset=mp_report.columns[1:])
        mp_report.to_csv(f'{out}/Metaproteomics/mp_entry_quant.tsv', sep='\t', index=False)
        shutil.copyfile(f'{out}/Metaproteomics/mp_entry_quant.tsv', f'{out}/Quantification/dea_input.tsv')


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')
    make_general_reports(snakemake.params.output, exps, did_assembly=snakemake.params.did_assembly)


if __name__ == '__main__':
    run()
