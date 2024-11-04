import pathlib
import pandas as pd
from time import gmtime, strftime
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

OUTPUT = config["output"]
EXPS = pd.DataFrame(config["experiments"])
pathlib.Path(f"{OUTPUT}").mkdir(parents=True, exist_ok=True)
EXPS.to_csv(f"{OUTPUT}/exps.tsv", sep = '\t', index = False)

mg_exps = EXPS[EXPS["Data type"] == 'dna']
mt_exps = EXPS[EXPS["Data type"] == 'mrna']
mp_exps = EXPS[EXPS["Data type"] == 'protein']

if len(mg_exps) == 0 and len(mt_exps) != 0:
    mg_exps = mt_exps
not_mp_exps = EXPS[EXPS["Data type"] != 'protein']
not_mg_exps = EXPS[EXPS["Data type"] != 'dna']
has_expression_data = len(not_mg_exps) > 0


def human_time(seconds):
    days = seconds // 86400
    if days > 0:
        return strftime(f"{days}d%Hh%Mm%Ss", gmtime(seconds))
    return strftime("%Hh%Mm%Ss", gmtime(seconds))


def sample_to_reads(wildcards):
    df = mg_exps[mg_exps['Sample'] == wildcards.sample].reset_index()
    return [f'{OUTPUT}/Preprocess/Trimmomatic/quality_trimmed_{df.iloc[row]["Name"]}{fr}.fq' for row in range(len(df))
           for fr in (['_forward_paired', '_reverse_paired'] if ',' in df.iloc[row]["Files"] else [''])]


def gene_calling_input(wildcards):
    df = mg_exps[mg_exps['Sample'] == wildcards.sample].reset_index()
    result = []
    for row in range(len(df)):
        if ',' in df.iloc[row]["Files"]:
            result.append(f'{OUTPUT}/Preprocess/Trimmomatic/{df.iloc[row]["Name"]}.assembled.fastq')
        else:
            result.append(f'{OUTPUT}/Preprocess/Trimmomatic/quality_trimmed_{df.iloc[row]["Name"]}.fq')
    return result