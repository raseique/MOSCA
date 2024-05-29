# -*- coding: utf-8 -*-
"""
MOSCA's Metaproteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

from glob import glob
import os
from pathlib import Path
import shutil
import pandas as pd
from tqdm import tqdm
import requests
from time import sleep
from mosca_tools import run_command, run_pipe_command, multiprocess_fun
import snakemake


class MetaproteomicsAnalyser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

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

    def add_reference_proteomes(self, upimapi_res, output, tax_lvl='genus'):
        taxids = pd.read_csv(upimapi_res, sep='\t', low_memory=False)[
            f'Taxonomic lineage IDs ({tax_lvl.upper()})'].dropna().unique().astype(int).tolist()
        with open(output, 'w') as f:
            for taxid in tqdm(taxids, desc=f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt'):
                try:
                    f.write(self.get_proteome_uniprot(taxid))
                except UnboundLocalError as e:
                    print(f'Failed to retrieve proteome for taxid {taxid}. Skipping.')

    def database_generation(
            self, mg_orfs, output, upimapi_res, contaminants_database=None, protease='Trypsin', threads=1,
            add_reference_proteomes=True, tax_lvl='genus'):
        """
        Build database from MG analysis
        :param mg_orfs:
        :param output:
        :param upimapi_res:
        :param contaminants_database:
        :param protease:
        :param threads:
        :param add_reference_proteomes:
        :param tax_lvl:
        """
        print(f'Generating new database in {output}')
        # Get reference proteomes for the various taxa
        if add_reference_proteomes:
            self.add_reference_proteomes(upimapi_res, f'{output}/ref_proteomes.fasta', tax_lvl=tax_lvl)
        # Add protease
        if protease == 'Trypsin':
            if not os.path.isfile(f'{output}/P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                run_command(f'wget https://www.uniprot.org/uniprot/P00761.fasta -P {output}')
            protease = f'{output}/P00761.fasta'
        else:  # is an inputed file
            if not os.path.isfile(protease):
                exit(f'Protease file does not exist: {protease}')
        files = [mg_orfs, protease]
        if add_reference_proteomes:
            files.append(f'{output}/ref_proteomes.fasta')
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
        for file in (
                ['database.fasta', 'predatabase.fasta'] + ['ref_proteomes.fasta'] if add_reference_proteomes else []):
            os.remove(f'{output}/{file}')

    def raw_to_mgf(self, file, out_dir, inside_container=False):
        """
        Convert raw file to mgf
        :param file:
        :param out_dir:
        :param inside_container:
        :return:
        """
        folder = os.path.dirname(os.path.abspath(file))
        file = os.path.basename(file)
        if inside_container:
            shutil.copyfile(f'{folder}/{file}', f'/data/{file}')
        run_pipe_command(
            f'docker run --rm -e WINEDEBUG=-all -v {"named_volume" if inside_container else folder}:/data '
            f'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/{file} --mgf '
            f'--filter "peakPicking cwt"')
        shutil.move(
            # if inside container, file is in the container's /data folder, else it will be in the original folder
            f'{"/data" if inside_container else folder}/{".".join(file.split(".")[:-1])}.mgf',
            f'{out_dir}/{".".join(file.split(".")[:-1])}.mgf')

    def spectra_in_proper_state(self, folder, out_dir, threads=1, inside_container=False):
        """
        Convert raw spectra files to MGF, and put all spectra in out_dir
        :param folder:
        :param out_dir:
        :param threads:
        :param inside_container:
        :return:
        """
        already_mgfs, files2convert = [], []
        for file in glob(f'{folder}/*'):
            if os.path.isfile(file):        # could be a folder
                if file.endswith('.mgf'):
                    already_mgfs.append(file)
                # elif file with termination replaced by mgf exists, skip
                elif os.path.isfile(f'{out_dir}/{".".join(file.split(".")[:-1])}.mgf'):
                    continue
                else:
                    files2convert.append(file)
        if len(files2convert) > 0:
            multiprocess_fun(
                self.raw_to_mgf, [(file, out_dir, inside_container) for file in files2convert], threads=threads)
        for file in set(already_mgfs):
            shutil.copyfile(f'{folder}/{file}', f'{out_dir}/{file}')

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
        Create decoy database
        :param database: str - Database to create decoy of
        """
        run_command(f'searchgui eu.isas.searchgui.cmd.FastaCLI -in {database} -decoy')

    def generate_parameters_file(
            self, output, protein_fdr=1, frag_tol=0.5, prec_tol=10, enzyme='Trypsin', mc=2,
            fixed_mods='Carbamidomethylation of C', variable_mods='Oxidation of M,Acetylation of protein N-term'):
        """
        param: output: name of parameters file
        param: protein_fdr: float - FDR at the protein level in percent
        param: frag_tol: float - fragment ion mass tolerance in ppm
        param: prec_tol: float - precursor ion mass tolerance in ppm
        param: enzyme: str - enzyme used for digestion
        param: mc: int - maximum number of missed cleavages
        param: fixed_mods: str - fixed modifications comma-separated
        param: variable_mods: str - variable modifications comma-separated
        returns: a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
        """
        fixed_mods, variable_mods = fixed_mods.split(','), variable_mods.split(',')
        run_pipe_command(
            f'''searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} -prec_tol {prec_tol} ''' +
            f'''-frag_tol {frag_tol} -enzyme {enzyme} -fixed_mods "{', '.join(fixed_mods)}" ''' +
            f'''-variable_mods "{', '.join(variable_mods)}" -mc {mc} -protein_fdr {protein_fdr}''')

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
            f"searchgui eu.isas.searchgui.cmd.SearchCLI -Xmx{max_memory}G -spectrum_files {spectra_folder} "
            f"-output_folder {output} -id_params {parameters_file} -threads {threads} -fasta_file {database} "
            f"{' '.join([f'-{engine} 1' for engine in search_engines])}")

    def browse_identification_results(
            self, spectra_folder, output, parameters_file, database, searchcli_output,  name, max_memory=4096):
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
                f'peptide-shaker -Xmx{max_memory}G eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files '
                f'{spectra_folder} -reference {name} -identification_files {searchcli_output} '
                f'-out {output} -fasta_file {database} -id_params {parameters_file}')
        except:
            print('Producing Peptide-Shaker result failed! Maybe no identifications were obtained?')

    def generate_reports(self, peptideshaker_output, reports_folder, reports=("10"), max_memory=40960):
        """
        Generates reports from output of PeptideShaker
        :param peptideshaker_output: peptideshaker output filename
        :param reports_folder: folder to where output reports
        :param reports: list of INTEGERS from 0 to 11 corresponding to the reports to output
        :param max_memory: maximum memory to use
        """
        print(f'Created {reports_folder}')
        Path(reports_folder).mkdir(parents=True, exist_ok=True)  # creates folder for reports
        run_command(
            f"peptide-shaker -Xmx{max_memory}G eu.isas.peptideshaker.cmd.ReportCLI -in {peptideshaker_output} "
            f"-out_reports {reports_folder} -reports {','.join(reports)}")

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

    def compomics_run(
            self, database, output, spectra_folders, name, params, threads=1, max_memory=4096, reports=('10')):
        """
        Run compomics workflow on the given spectra folders
        :param params:
        :param database:
        :param output:
        :param spectra_folders:
        :param name:
        :param threads:
        :param max_memory:
        :param reports:
        :return:
        """
        self.peptide_spectrum_matching(
            spectra_folders, output, params, database, threads=threads, max_memory=max_memory)
        self.browse_identification_results(
            spectra_folders, f'{output}/ps_output.psdb', params, database, f'{output}/searchgui_out.zip',
            name, max_memory=max_memory)
        try:
            self.generate_reports(
                f'{output}/ps_output.psdb', f'{output}/reports', max_memory=max_memory,
                reports=reports)
        except Exception as e:
            print(f'Reporter generation sent an error: {e}\nBut likely worked!')

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

    def run(self):
        Path(snakemake.params.output).mkdir(parents=True, exist_ok=True)
        # 1st database construction
        self.database_generation(
            snakemake.params.mg_db, snakemake.params.output, snakemake.params.up_res,
            contaminants_database=snakemake.params.contaminants_database,
            protease=snakemake.params.protease, add_reference_proteomes=snakemake.params.add_reference_proteomes,
            tax_lvl=snakemake.params.taxa_lvl)
        self.create_decoy_database(f'{snakemake.params.output}/1st_search_database.fasta')
        self.split_database(
            f'{snakemake.params.output}/1st_search_database_concatenated_target_decoy.fasta', n_proteins=5000000)
        try:  # try/except - https://github.com/compomics/searchgui/issues/217
            self.generate_parameters_file(f'{snakemake.params.output}/1st_params.par', protein_fdr=100)
        except:
            print('An illegal reflective access operation has occurred. But MOSCA can handle it.')

        # 2nd database construction - peak-pick spectra and get proteins from taxa with at least one identification
        up_res = pd.read_csv(f'real_mgmp_output/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t', low_memory=False)

        proteins_identified = []
        tax_lvl = snakemake.params.tax_level
        for i in range(len(snakemake.params.names)):
            out = f'{snakemake.params.output}/{snakemake.params.names[i]}'
            for foldername in ['spectra', '2nd_search']:
                Path(f'{out}/{foldername}').mkdir(parents=True, exist_ok=True)
            self.spectra_in_proper_state(
                snakemake.params.spectra_folders[i], f'{out}/spectra/',
                inside_container=snakemake.params.inside_container)
            j = 1
            for database in sorted(glob(f'{snakemake.params.output}/1st_search_database_*.part_*.fasta')):
                Path(f'{out}_{j}/1st_search').mkdir(parents=True, exist_ok=True)
                self.compomics_run(
                    database, f'{out}_{j}/1st_search', f'{out}/spectra', snakemake.params.names[i],
                    f'{snakemake.params.output}/1st_params.par', threads=snakemake.threads,
                    max_memory=snakemake.params.max_memory, reports=('4'))
                df = pd.read_csv(
                    f'{out}_{j}/1st_search/reports/{snakemake.params.names[i]}_'
                    f'Default_PSM_Report_with_non-validated_matches.txt',
                    sep='\t', index_col=0)
                for protein_group in df['Protein(s)'].str.split(','):
                    proteins_identified += protein_group
                j += 1
        identified_taxa = up_res[up_res[f'Taxonomic lineage ({tax_lvl})'].notnull() & up_res['qseqid'].isin(
            proteins_identified)][f'Taxonomic lineage ({tax_lvl})'].unique()
        with open(f'{snakemake.params.output}/2nd_search_proteins.txt', 'w') as f:
            # write proteins corresponding to identified taxa
            f.write('\n'.join(set(up_res[up_res[f'Taxonomic lineage ({tax_lvl})'].isin(identified_taxa)]['qseqid'])))
        run_command(
            f'seqkit grep -f {snakemake.params.output}/2nd_search_proteins.txt -o '
            f'{snakemake.params.output}/2nd_search_database.fasta {snakemake.params.output}/1st_search_database.fasta')

        self.create_decoy_database(f'{snakemake.params.output}/2nd_search_database.fasta')
        self.generate_parameters_file(f'{snakemake.params.output}/2nd_params.par', protein_fdr=1)

        # TODO - check if splitting the database will be necessary
        #self.split_database(f'{args.output}/2nd_search_database_concatenated_target_decoy.fasta', n_proteins=5000000)

        # Protein identification and quantification
        spectracounts = pd.DataFrame(columns=['Main Accession'])
        for name in snakemake.params.names:
            out = f'{snakemake.params.output}/{name}'
            self.compomics_run(
                f'{snakemake.params.output}/2nd_search_database_concatenated_target_decoy.fasta',
                f'{out}/2nd_search', f'{out}/spectra', name, f'{snakemake.params.output}/2nd_params.par',
                threads=snakemake.threads, max_memory=snakemake.params.max_memory, reports=['9', '10'])
            spectracounts = pd.merge(spectracounts, pd.read_csv(
                f'{out}/2nd_search/reports/{name}_Default_Protein_Report.txt', sep='\t', index_col=0
            )[['Main Accession', '#PSMs']].rename(columns={'#PSMs': name}), how='outer', on='Main Accession')
        spectracounts.groupby('Main Accession')[snakemake.params.names].sum().reset_index().fillna(value=0.0).to_csv(
            f'{snakemake.params.output}_mp.spectracounts', sep='\t', index=False)


if __name__ == '__main__':
    MetaproteomicsAnalyser().run()
