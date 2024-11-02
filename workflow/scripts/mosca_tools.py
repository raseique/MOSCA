# -*- coding: utf-8 -*-
"""
General tools for support of MOSCA's functionalities

By João Sequeira

Jun 2017
"""

from glob import glob
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
from subprocess import run, Popen, PIPE, check_output
import sys
import time

blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
              'evalue', 'bitscore']


def run_command(bashCommand, output='', mode='w', sep=' ', print_message=True, verbose=True):
    if print_message:
        print(f"{bashCommand.replace(sep, ' ')}{' > ' + output if output != '' else ''}")
    if output == '':
        run(bashCommand.split(sep), stdout=sys.stdout if verbose else None, check=True)
    else:
        with open(output, mode) as output_file:
            run(bashCommand.split(sep), stdout=output_file)


def run_pipe_command(bashCommand, output='', mode='w', sep=' ', print_message=True):
    if print_message:
        print(bashCommand)
    if output == '':
        Popen(bashCommand, stdin=PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return Popen(bashCommand, stdin=PIPE, shell=True, stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            Popen(bashCommand, stdin=PIPE, shell=True, stdout=output_file).communicate()


def parse_blast(blast):
    result = pd.read_csv(blast, sep='\t', header=None)
    result.columns = blast_cols
    return result


def check_bowtie2_index(index_prefix):
    files = glob(index_prefix + '*.bt2')
    if len(files) < 6:
        return False
    return True


def generate_mg_index(reference, index_prefix):
    run_pipe_command(f'bowtie2-build {reference} {index_prefix} 1> {index_prefix}.log 2> {index_prefix}.err')


def align_reads(reads, index_prefix, sam, report, log=None, threads=6):
    run_pipe_command(
        f'bowtie2 -x {index_prefix} -1 {reads[0]} -2 {reads[1]} -S {sam} -p {threads} 1> {report} 2> {log}')


def parse_fasta(file):
    print(f'Parsing {file}')
    lines = [line.rstrip('\n') for line in open(file)]
    i = 0
    sequences = dict()
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:]
            sequences[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                sequences[name] += lines[i]
                i += 1
    return sequences


def build_gff_from_contigs(contigs, output, assembler=None):
    gff = pd.DataFrame()
    contigs = parse_fasta(contigs)
    contigs = [[k, v] for k, v in contigs.items()]
    gff["seqid"] = [contig[0] for contig in contigs]
    size = len(contigs)
    gff["source"] = ['.' if assembler is None else assembler] * size
    gff["type"] = ['exon'] * size
    gff["start"] = ['1'] * size
    gff["end"] = [len(contigs[1]) for contigs in contigs]
    gff["score"] = ['.'] * size
    gff["strand"] = ['.'] * size
    gff["phase"] = ['.'] * size
    gff["attributes"] = ['gene_id=' + contig[0] for contig in contigs]
    gff.to_csv(output, sep='\t', index=False, header=False)


def build_gff_from_orfs(orfs, output):
    gff = pd.DataFrame()
    qseqids = parse_fasta(orfs).keys()
    parts = [qid.split('_') for qid in qseqids]
    gff["seqid"] = qseqids
    size = gff.size
    gff["source"] = ['UniProtKB'] * size
    gff["type"] = ['exon'] * size
    gff["start"] = [1] * size
    gff["end"] = [int(part[-2]) - int(part[-3]) + 1 for part in parts]
    gff["score"] = ['.'] * size
    gff["strand"] = [part[-1] for part in parts]
    gff["phase"] = ['.'] * size
    gff["Name"] = [f"Name={'_'.join(part)}" for part in parts]
    gff.to_csv(output, sep='\t', index=False, header=False)


def perform_alignment(reference, reads, basename, threads=1):
    """
    Input:
        reference: name of contigs file from assembly
        reads: list, [forward reads, reverse reads]
        basename: basename of outputs
        threads: number of threads to use
        attribute: str - identifier to group quantification for
        type_of_reference: str - 'contigs' or 'orfs'
    Output:
        Will generate a bowtie2 index named contigs.replace(.fasta,_index),
        GFF annotation file named blast.replace(.blast,.gff)
        SAM alignment and READCOUNTS files named output + .sam and .readcounts
    """
    ext = f".{reference.split('.')[-1]}"
    if not check_bowtie2_index(reference.replace(ext, '_index')):
        print('INDEX files not found. Generating new ones')
        generate_mg_index(reference, reference.replace(ext, '_index'))
    else:
        print(f"INDEX was located at {reference.replace(ext, '_index')}")
    if not os.path.isfile(f'{basename}.log'):
        align_reads(reads, reference.replace(ext, '_index'), f'{basename}.sam',
                    f'{basename}_bowtie2_report.txt', log=f'{basename}.log', threads=threads)
    else:
        print(f'{basename}.log was found!')
    run_pipe_command(
        f"""samtools view -F 260 -S {basename}.sam | cut -f 3 | sort | uniq -c | \
        awk '{{printf("%s\\t%s\\n", $2, $1)}}'""", output=f'{basename}.readcounts')


def fastq2fasta(fastq, output):
    run_command(f"""paste - - - - < {fastq} | cut -f 1,2 | sed 's/^@/>/' | tr \"\t\" \"\n\" > {output}""")


def timed_message(message):
    print(f'{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}: {message}')


def normalize_counts_by_size(readcounts, reference):
    run_pipe_command(
        f"seqkit fx2tab {reference} | sort -k 1b,1 | awk '{{print $1\"\\t\"length($2)}}' | join - {readcounts} | "
        f"awk '{{print $1\"\\t\"$3/$2}}'", output=readcounts.replace('.readcounts', '.readcounts.norm'))


def add_abundance(data, readcounts, name, origin_of_data='metagenomics'):
    """
    Input:
        data: pd.DataFrame - Protein Report on the making
        readcounts: readcounts file from quantification with htseq-count
        origin_of_data: 'metagenomics', 'metatranscriptomics' defines column of joining
        name: name of column to add to relation
    Output:
        'data' df will receive additional column with abundance/expression
        information. This column will be named 'name'
    """
    readcounts = pd.read_csv(readcounts, sep='\t', header=None, names=['qseqid', name])
    readcounts[name] = readcounts[name].fillna(value=0)
    if origin_of_data == 'metagenomics':
        readcounts['Contig'] = [qseqid.split('_')[1] for qseqid in readcounts['qseqid']]
        del readcounts['qseqid']
        return pd.merge(data, readcounts, on='Contig', how='left')
    elif origin_of_data == 'metatranscriptomics':
        return pd.merge(data, readcounts, on='qseqid', how='left')
    elif origin_of_data == 'compomics':
        pass


def multi_sheet_excel(output, data, sheet_name='Sheet', max_lines=1000000, index=False):
    k = 1
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    if len(data) < max_lines:
        data.to_excel(writer, sheet_name=f'{sheet_name}', index=index)
    else:
        for i in range(0, len(data), max_lines):
            j = min(i + max_lines, len(data))
            data.iloc[i:(i + j)].to_excel(writer, sheet_name=f'{sheet_name} ({k})', index=index)
            k += 1
    writer.close()


def expand_by_list_column(df, column='Pathway'):
    lens = [len(item) for item in df[column]]
    dictionary = dict()
    for column in df.columns:
        dictionary[column] = np.repeat(df[column].values, lens)
    dictionary[column] = np.concatenate(df[column].values)
    return pd.DataFrame(dictionary)


def parse_fastqc_report(filename):
    data = {}
    file = open(filename).read().split('\n')
    i = 1
    while i < len(file):
        if file[i].startswith('>>') and file[i] != '>>END_MODULE':
            name, flag = file[i][2:].split('\t')[0], file[i][2:].split('\t')[1]
            if name == 'Sequence Duplication Levels':
                i += 1
            i += 1
            if file[i] == '>>END_MODULE':
                data[name] = (flag, pd.DataFrame())
            else:
                labels = file[i][1:].split('\t')
                i += 1
                partial_data = np.array(labels)
                while i < len(file) and not file[i].startswith('>>'):
                    partial_data = np.append(partial_data, file[i].split('\t'))
                    i += 1
                partial_data = np.reshape(partial_data, (int(partial_data.size / len(labels)), len(labels)))
                data[name] = (flag, pd.DataFrame(data=partial_data[1:, 1:], index=partial_data[1:, 0],
                                                 columns=partial_data[0, 1:]))
        i += 1
    return data


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def count_lines(file):
    return int(check_output(f"wc -l {file}", shell=True).decode('utf8').split()[0])


def sort_alphanumeric(alphanumeric_list):
    return sorted(alphanumeric_list, key=lambda item: (
        int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))


def generate_expression_matrix(readcount_files, header, output):
    expression_matrix = pd.DataFrame()
    for file in readcount_files:
        df = pd.read_csv(file, sep='\t', index_col=0, header=None)
        df.columns = [file.split('/')[-1].rstrip('.readcounts')]
        expression_matrix = pd.merge(expression_matrix, df, how='outer',
                                     left_index=True, right_index=True)
    expression_matrix = expression_matrix[1:-5]  # remove non identified proteins (*) and the metrics at the end
    expression_matrix = expression_matrix.fillna(value=0).astype(
        int)  # set not identified proteins expression to 0, and convert floats, since DeSEQ2 only accepts integers
    expression_matrix.index.name = 'geneid'
    expression_matrix.columns = header
    expression_matrix.to_csv(output, sep='\t')


def fastqc_name(filename):
    return filename.replace("stdin:", "").replace(".gz", "").replace(".bz2", "").replace(".txt", "").replace(
        ".fastq", "").replace(".fq", "").replace(".csfastq", "").replace(".sam", "").replace(".bam", "")


def multiprocess_fun(fun, args, threads):
    """
    Run a function in parallel using multiprocessing
    :param fun:
    :param args: list of tuples
    :param threads:
    :return:
    """
    with Pool(processes=threads) as p:
        p.starmap(fun, args)