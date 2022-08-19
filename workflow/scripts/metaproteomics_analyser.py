# -*- coding: utf-8 -*-
"""
MOSCA's Metaproteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

from glob import glob
import os
import pathlib
import shutil

import pandas as pd
import argparse
from lxml import etree
from mosca_tools import run_command, sort_alphanumeric, parse_blast, run_pipe_command, multiprocess_fun, count_on_file
from tqdm import tqdm
import requests
from time import sleep


class MetaproteomicsAnalyser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's metaproteomics analysis")
        parser.add_argument("-sfs", "--spectra-folders", help="List of folders with spectra to be analysed")
        parser.add_argument("-ns", "--names", help="List of names for datasets")
        parser.add_argument("-t", "--threads", help="Number of threads to use.")
        parser.add_argument("-o", "--output", help="Output directory, where results are stored.")
        parser.add_argument("-w", "--workflow", help="Workflow to use", choices=['maxquant', 'compomics'])
        parser.add_argument(
            "-db", "--database", help="Database file (FASTA format) from metagenomics for protein identification")
        parser.add_argument(
            "-cdb", "--contaminants-database", default=None,
            help="Database file (FASTA format) with contaminant sequences")
        parser.add_argument("-mpar", "--metaphlan-result", help="Results from MetaPhlan taxonomic annotation")
        parser.add_argument(
            "-rtl", "--references-taxa-level", default='genus',
            choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
            help="Taxonomic level to retrieve reference proteomes from")
        parser.add_argument(
            "-bs", "--batch-size", type=int, default=5000, help="How many IDs to submit per request to NCBI")
        parser.add_argument("-ma", "--max-attempts", default=3, help="Maximum attempts to access NCBI")
        parser.add_argument("--protease", help="Filename in fasta format of protease sequence", default='Trypsin')
        parser.add_argument(
            "-mmem", "--max-memory", type=int, default=4, help="Maximum memory (Gb) to use for Peptide-to-Spectrum matching [4]")
        parser.add_argument(
            "-rd", "--resources-directory", default=os.path.expanduser('~/resources'),
            help="Directory for storing databases and other important files")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        args.spectra_folders = args.spectra_folders.split(',')
        args.names = args.names.split(',')
        if len(args.names) != len(args.spectra_folders):
            exit('Length of names and spectra_folders must be equal.')
        args.max_memory *= 1024
        return args

    def get_proteome_uniprot(self, taxid, max_tries=3):
        tries = 0
        done = False
        while tries < max_tries and not done:
            try:
                res = requests.get(
                    f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A{taxid}%29')
                done = True
            except:
                print(f'Failed! {max_tries - tries} tries remaining.')
                tries += 1
                sleep(10)
        return res.content.decode('utf8')

    def add_reference_proteomes(self, mpa_result, output, references_taxa_level='genus'):
        """
        Input:
            mpa_result: str - filename of MetaPhlan result
            output: str - name of folder to output reference proteomes
            references_taxa_level: str - one of 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
            max_attemps: int - number of requests failed before giving up
        """
        taxa_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        mpa_data = pd.read_csv(mpa_result, sep='\t', skiprows=3)
        mpa_data = mpa_data[mpa_data['NCBI_tax_id'].str.count('\|') == taxa_levels.index(references_taxa_level)]
        taxids = [ide.split('|')[-1] for ide in mpa_data['NCBI_tax_id']]
        with open(output, 'w') as f:
            for taxid in tqdm(taxids, desc=f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt'):
                f.write(self.get_proteome_uniprot(taxid))

    def database_generation(
            self, mg_orfs, output, mpa_result, contaminants_database=None, protease='Trypsin',
            references_taxa_level='genus', threads=1):
        """
        Build database from MG analysis
        :param mg_orfs:
        :param output:
        :param mpa_result:
        :param contaminants_database:
        :param protease:
        :param references_taxa_level:
        :param threads:
        """
        print(f'Generating new database in {output}')
        # Get reference proteomes for the various taxa
        self.add_reference_proteomes(
            mpa_result, f'{output}/ref_proteomes.fasta', references_taxa_level=references_taxa_level)
        # Add protease
        if protease == 'Trypsin':
            if not os.path.isfile(f'{output}/P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                run_command(f'wget https://www.uniprot.org/uniprot/P00761.fasta -P {output}')
            protease = f'{output}/P00761.fasta'
        else:  # is an inputed file
            if not os.path.isfile(protease):
                exit(f'Protease file does not exist: {protease}')
        files = [mg_orfs, f'{output}/ref_proteomes.fasta', protease]
        if contaminants_database is not None:
            self.verify_crap_db(contaminants_database)
            files.append(contaminants_database)
        run_command(f"cat {' '.join(files)}", output=f'{output}/predatabase.fasta', mode='w')
        # Join aminoacid lines, and remove empty lines
        run_pipe_command(f"awk '{{if ($0 ~ /^>/) {{print \"\\n\"$0}} else {{printf $0}}}}' {output}/predatabase.fasta",
                         output=f'{output}/database.fasta')
        run_pipe_command(f"awk '{{print $1}}' {output}/database.fasta", output=f'{output}/predatabase.fasta')
        # Remove asterisks (non identified aminoacids) and plicas
        run_pipe_command(f"""sed -i "s/[*\']//g" {output}/predatabase.fasta""")
        # Remove duplicate sequences
        run_command(
            f'seqkit rmdup -s -i -w 0 -o {output}/1st_search_database.fasta -D {output}/seqkit_duplicated.detail.txt '
            f'-j {threads} {output}/predatabase.fasta')
        for file in ['database.fasta', 'predatabase.fasta', 'unique.fasta', 'ref_proteomes.fasta']:
            os.remove(f'{output}/{file}')

    def raw_to_mgf(self, file, out_dir):
        """
        Convert raw file to mgf
        :param file:
        :param out_dir:
        :param outfmt:
        :param peak_picking:
        :return:
        """
        folder, filename = os.path.split(file)
        pathlib.Path('named_volume').mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f'{folder}/{filename}', f'/data/{filename}')
        run_pipe_command(
            f'docker run --rm -e WINEDEBUG=-all -v named_volume:/data '
            f'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/{filename} --mgf '
            f'--filter "peakPicking cwt"')
        shutil.copyfile(f'/data/{".".join(filename.split(".")[:-1])}.mgf', out_dir)

    def spectra_in_proper_state(self, folder, out_dir, threads=1):
        """
        Convert raw spectra files to MGF, and put all spectra in out_dir
        :param threads:
        :param folder:
        :param out_dir:
        :return:
        """
        already_mgfs, files2convert = [], []
        for file in glob(f'{folder}/*'):
            if os.path.isfile(file):
                if file.endswith('.mgf'):
                    already_mgfs.append(file)
                else:
                    files2convert.append(file)
        if len(files2convert) > 0:
            multiprocess_fun(self.raw_to_mgf, [(file, out_dir) for file in files2convert], threads=threads)
        for file in already_mgfs:
            os.copy(f'{folder}/{file}', f'{out_dir}/{file}')

    def verify_crap_db(self, contaminants_database='MOSCA/Databases/metaproteomics/crap.fasta'):
        """
        Checks if contaminants database exists. If not, cRAP will be downloaded.
        :param contaminants_database:
        :return:
        """
        if os.path.isfile(contaminants_database):
            print(f'cRAP database exists at {contaminants_database}')
        else:
            print(f'cRAP database not found at {contaminants_database}. Downloading cRAP database.')
            run_command(f'wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta -O {contaminants_database}')

    def create_decoy_database(self, database):
        """
        input:
            protein_fasta: fasta file with proteins from MG analysis, plus trypsin
            and cRAP sequences
        output:
            a FASTA file named "database + _concatenated_target_decoy.fasta" will be
            created, containing interleaved original and decoy sequences
        """
        decoy_database = database.replace('.fasta', '_concatenated_target_decoy.fasta')
        if not os.path.isfile(decoy_database):
            run_command(f'searchgui eu.isas.searchgui.cmd.FastaCLI -in {database} -decoy')
        else:
            print(f'{decoy_database} already exists!')

    def generate_parameters_file(self, output, protein_fdr=1):
        """
        input:
            output: name of parameters file
            database: name of FASTA decoy database
            protein_fdr: float - FDR at the protein level in percent
        output:
            a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
        """
        run_pipe_command(
            f'searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} -prec_tol 10 '
            f'-frag_tol 0.02 -enzyme Trypsin -fixed_mods "Carbamidomethylation of C" -variable_mods '
            f'"Oxidation of M, Acetylation of protein N-term" -mc 2 -protein_fdr {protein_fdr}')

    def split_database(self, database, n_proteins=5000000):
        """
        Split database into smaller files
        :param database:
        :param n_proteins:
        :return:
        """
        run_command(f'seqkit split -s {n_proteins} {database} -O {os.path.dirname(database)}')

    def peptide_spectrum_matching(
            self, spectra_folder, output, parameters_file, database, threads=12, max_memory=4096,
            search_engines=('xtandem', 'myrimatch', 'msgf')):
        """
        input:
            spectra_folder: folder containing the raw spectra files
            output: folder to output results
            parameters_file: parameters filename
            search_engines: search engines to perform PSM
        output:
            a "searchgui_out.zip" file will be created in the output folder
        """
        run_command(
            f"searchgui eu.isas.searchgui.cmd.SearchCLI -Xmx{max_memory}M -spectrum_files {spectra_folder} "
            f"-output_folder {output} -id_params {parameters_file} -threads {threads} -fasta_file {database} "
            f"{' '.join([f'-{engine} 1' for engine in search_engines])}")

    def browse_identification_results(
            self, spectra_folder, fasta_file, searchcli_output, peptideshaker_output, name, max_memory=4096):
        """
        input:
            spectra_folder: folder containing the raw spectra files
            output: folder to output results
            parameters_file: parameters filename
            searchcli_output: searchcli output filename
            peptideshaker_output: peptideshaker output filename
            experiment_name: name of experiment
            sample_name: name of sample
            replicate_number: number of replicate (STRING)
        output:
            a file will be outputed with the validation of the PSMs and protein
            identifications
        """
        try:
            run_command(
                f'peptide-shaker -Xmx{max_memory}M eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files '
                f'{spectra_folder} -reference {name} '
                f'-identification_files {searchcli_output} -out {peptideshaker_output} -fasta_file {fasta_file}')
        except:
            print('Producing Peptide-Shaker result failed! Maybe no identifications were obtained?')

    def generate_reports(self, peptideshaker_output, reports_folder, reports_list=[10]):
        """
        input:
            peptideshaker_output: peptideshaker output filename
            reports_folder: folder to where output reports
            reports_list: list of INTEGERS from 0 to 11 corresponding to the reports
            to output
        output:
            if it doesn't exist, "reports_folder" will be created
            reports will be outputed to "reports_folder"
        """
        print(f'Created {reports_folder}')
        pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)  # creates folder for reports
        run_command(
            f"peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {peptideshaker_output} "
            f"-out_reports {reports_folder} -reports {','.join(reports_list)}")

    def join_ps_reports(self, files, local_fdr=None, validation=False):
        result = pd.DataFrame(columns=['Main Accession'])
        for file in files:
            data = pd.read_csv(file, sep='\t', index_col=0)
            if local_fdr is not None:
                data = data[data['Confidence [%]'] > 100 - local_fdr]
            if validation:
                data = data[data['Validation'] == 'Confident']
            if len(data) > 0:
                name = file.split('/')[-1].split('_')[1]
                data = data[['Main Accession', '#PSMs']]
                data.columns = ['Main Accession', name]
                result = pd.merge(result, data, on='Main Accession', how='outer')
        return result

    def spectra_counting(self, protein_reports, output, blast=None, uniprot_ids=False, samples_names=None):
        """
        input:
            protein_report: name of file containing protein report from PeptideShaker
            output: name of file to output
            blast: name of blast file with UniProt annotations - only needed if uniprot_ids = False
            uniprot_ids: boolean, True if protein sequences already are named with UniProt IDs
                If False, will need blast file to retrieve the UniProt IDs from
            samples_names: names of samples respective to protein_reports
        output:
            A tab separated spectra count file named "output", with UniProt IDs
            and corresponding spectra quantification
        """
        protein_reports = sort_alphanumeric(protein_reports)
        if samples_names is None:
            samples_names = [filename.split('/')[-3] for filename in
                             protein_reports]  # samples names are the folder containing the reports folder
        spectra_count = pd.DataFrame(columns=['Main Accession'])
        for i in range(len(protein_reports)):
            report = pd.read_csv(protein_reports[i], sep='\t', index_col=0)
            if not uniprot_ids:
                blast = parse_blast(blast)
                blast['sseqid'] = [ide.split('|')[-1] for ide in blast.sseqid]
                report = pd.merge(report, blast[['qseqid', 'sseqid']], left_on='Main Accession', right_on='qseqid')
                report = report[['sseqid', '#PSMs']]
            else:
                report = report[['Main Accession', '#PSMs']]
            report.columns = ['Main Accession', samples_names[i]]
            report = report.groupby('Main Accession')[samples_names[i]].sum().reset_index()
            spectra_count = pd.merge(spectra_count, report, on='Main Accession', how='outer')
        spectra_count[samples_names] = spectra_count[samples_names].fillna(value=0).astype(int)
        spectra_count.to_csv(output, sep='\t', index=False)

    def compomics_workflow(self, mg_db, output, spectra_folders, name, threads=1, protein_fdr=1, max_memory=4096):
        """
        Run compomics workflow on the given spectra folders
        :param mg_db:
        :param output:
        :param spectra_folders:
        :param name:
        :param threads:
        :param protein_fdr:
        :param max_memory:
        :return:
        """
        """
        self.create_decoy_database(mg_db)
        try:  # try/except - https://github.com/compomics/searchgui/issues/217
            self.generate_parameters_file(f'{output}/params.par', protein_fdr=protein_fdr)
        except:
            print('An illegal reflective access operation has occurred. But MOSCA can handle it.')
        """
        self.split_database(mg_db.replace('.fasta', '_concatenated_target_decoy.fasta'), n_proteins=5000000)
        for database in glob(f'{os.path.dirname(mg_db)}.part_*.fasta'):
            self.peptide_spectrum_matching(
                spectra_folders, output, f'{output}/params.par', database, threads=threads, max_memory=max_memory)
            '''
            self.browse_identification_results(
                spectra_folders, f'{output}/params.par', f'{output}/searchgui_out.zip', f'{output}/ps_output.cpsx', 
                name, max_memory=max_memory)
            try:  # try/except - if no identifications are present, will throw an error
                self.generate_reports(f'{output}/ps_output.cpsx', f'{output}/reports')
            except:
                print('No identifications?')
            self.spectra_counting(
                f'{output}/reports/{experiment_name}_{sample_name}_{replicate_number}_Default_Protein_Report.txt',
                self.blast, f'{self.output}/Spectra_counting.tsv')
            '''

    """
    Metaproteomics with MaxQuant
    """

    def create_mqpar(self, output):
        """
        input:
            output: name of standard parameters file to create
        output:
            a standard parameters file for MaxQuant named "output" will be created
        """
        if os.path.isfile(output):  # the create file command will not create a new one if the file already exists
            os.remove(
                output)  # even if that file already has non-default information in it, messing with the next commands
        run_command(f'maxquant {output} --create')

    def edit_maxquant_mqpar(
            self, mqpar, mg_db, spectra_folders, experiment_names, threads=1, spectra_format='RAW', protein_fdr=0.01):
        """
        input:
            mqpar: name of the mqpar.xml file to have its values changed
            fasta_database: name of the FASTA database for metaproteomics
            spectra_folder: name of the folder with RAW spectra files
            experiment_names: list of experiments, with one element for each file
        output:
            the "file" file will be updated with the new parameters
        """
        print('Updating parameters file information.')
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(mqpar, parser)
        root = tree.getroot()
        root.find("fastaFiles/FastaFileInfo/fastaFilePath").text = mg_db
        print(f'Fasta database = {mg_db}')
        root.find("separateLfq").text = 'True'
        root.find("numThreads").text = str(threads)
        root.find("proteinFdr").text = str(protein_fdr)
        print(f'Number of threads = {threads}')
        for child in [
            'filePaths/string', 'experiments/string', 'fractions/short', 'ptms/boolean', 'paramGroupIndices/int']:
            tree.xpath(child)[0].getparent().remove(tree.xpath(child)[0])  # removes the old information
        filePaths = root.find("filePaths")
        experiments = root.find("experiments")
        fractions = root.find("fractions")
        ptms = root.find("ptms")
        paramGroupIndices = root.find("paramGroupIndices")
        for i in range(len(spectra_folders)):
            files = sort_alphanumeric(glob(f'{spectra_folders[i]}/*.{spectra_format}'))
            print(f'Adding files from folder: {spectra_folders[i]}')
            for file in files:
                etree.SubElement(filePaths, 'string').text = file
                etree.SubElement(experiments, 'string').text = experiment_names[i]
                etree.SubElement(fractions, 'short').text = '32767'
                etree.SubElement(ptms, 'boolean').text = 'True'
                etree.SubElement(paramGroupIndices, 'int').text = '0'
        root.find("parameterGroups/parameterGroup/lfqMode").text = '1'
        tree.write(mqpar, pretty_print=True)
        print(f'Parameters file is available at {mqpar}')

    def run_maxquant(self, mqpar, spectra_folder, output_folder):
        """
        Check if MaxQuant's output folder exist and delete it, run MaxQuant, and move output folder to desired place
        :param mqpar:
        :param spectra_folder:
        :param output_folder:
        """
        for directory in [f'{spectra_folder}/combined', output_folder]:
            if os.path.isdir(directory):
                shutil.rmtree(directory, ignore_errors=True)
        run_command(f'maxquant {mqpar}')
        os.rename(f'{spectra_folder}/combined', output_folder)

    def maxquant_workflow(
            self, mqpar, mg_db, spectra_folders, experiment_names, output, threads=1, spectra_format='RAW',
            protein_fdr=0.01):
        """
        Input:
            mqpar: str - filename of parameters file to create
            database: str - filename of protein database
            spectra_folder: str - name of folder containing spectra
            experiment_names: list - of names to use as experiment names
            output: str - name of folder to output maxquant results
            threads: int - number of threads to use
        """
        self.create_mqpar(mqpar)
        self.edit_maxquant_mqpar(
            mqpar, mg_db, spectra_folders, experiment_names, threads=threads, spectra_format=spectra_format,
            protein_fdr=protein_fdr)
        self.run_maxquant(mqpar, spectra_folders, output)

    def select_proteins_for_second_search(self, original_db, output, results_files, column='Main Accession'):
        proteins = []
        for file in results_files:
            proteins += pd.read_csv(file, sep='\t')[column].tolist()
        proteins = set(proteins)
        print(f'Selected {len(proteins)} proteins for 2nd peptide-to-spectrum matching.')
        with open(f'{output}/2nd_search_ids.txt', 'w') as f:
            f.write('\n'.join(proteins))
        # the '\|(.*)\|' rule also works for the selection of full ids (like the ones coming from FragGeneScan)
        run_pipe_command(
            f"seqkit grep {original_db} -w 0 --id-regexp '\|(.*)\|' -f {output}/2nd_search_ids.txt",
            output=f'{output}/2nd_search_database.fasta')

    def background_inputation(self, df):
        return df.fillna(value=df.min().min())

    def censored_inputation(self, df, replicates):
        for replicate in replicates:
            for line in df.index:
                if df.loc[df.index[0]][replicates[0]].isnull().sum() > 1:
                    df.loc[df.index[0]][replicates[0]] = (
                        self.background_inputation(df.loc[df.index[0]][replicates[0]]))

    def run(self):
        args = self.get_arguments()
        '''
        self.database_generation(
            args.database, args.output, args.metaphlan_result,
            contaminants_database=args.contaminants_database, protease=args.protease,
            references_taxa_level=args.references_taxa_level, max_attempts=args.max_attemps)
        '''
        for i in range(len(args.names)):
            out = f'{args.output}/{args.names[i]}'
            print(out)
            for foldername in ['spectra', '1st_search', '2nd_search']:
                pathlib.Path(f'{out}/{foldername}').mkdir(parents=True, exist_ok=True)

            if args.workflow == 'maxquant':
                self.maxquant_workflow(
                    f'{args.output}/mqpar.xml', f'{args.output}/1st_search_database.fasta',
                    args.spectra_folders.split(','), args.experiment_names.split(','), args.output,
                    threads=args.threads, spectra_format='RAW', protein_fdr=1)

            elif args.workflow == 'compomics':
                '''
                self.spectra_in_proper_state(args.spectra_folders[i], f'{out}/spectra')
                '''
                self.compomics_workflow(
                    f'{args.output}/1st_search_database.fasta', f'{out}/1st_search', f'{out}/spectra', args.names[i],
                    threads=args.threads, protein_fdr=100, max_memory=args.max_memory)
            '''
            self.select_proteins_for_second_search(
                f'{args.output}/Metaproteomics/{sample}/1st_search_database.fasta',
                f'{args.output}/Metaproteomics/{sample}',
                glob(f'{args.output}/Metaproteomics/{sample}/1st_search/*/reports/'
                     '*_Protein_Report_with_non-validated_matches.txt'), column='Main Accession')
            for name in set(partial_sample['Name']):
                self.compomics_workflow(
                    f'{args.output}/Metaproteomics/{sample}/2nd_search_database.fasta',
                    f'{args.output}/Metaproteomics/{sample}/2nd_search/{name}',
                    f'{args.output}/Metaproteomics/{sample}/spectra/{name}',
                    threads=args.threads, protein_fdr=1, max_memory=args.max_memory,
                    experiment_name=sample, sample_name=name)
            self.join_ps_reports(
                ['{0}/Metaproteomics/{1}/2nd_search/{2}/reports/{1}_{2}_1_Default_Protein_Report.txt'.format(
                    args.output, sample, name) for name in set(partial_sample['Name'])],
                local_fdr=5, validation=False
            ).to_csv(f'{args.output}/Metaproteomics/{sample}/quantification_table.tsv', sep='\t')
            '''
        else:
            print('Not a valid workflow option!')


if __name__ == '__main__':
    MetaproteomicsAnalyser().run()
