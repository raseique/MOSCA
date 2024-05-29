rule binning:
    input:
        reads = expand("{output}/Preprocess/{{sample}}{fr}.fastq", output=OUTPUT,
            fr=(['_forward', '_reverse'] if EXPS["Files"].str.contains(',').tolist() else '')),
        contigs = expand("{output}/Assembly/{{sample}}/scaffolds.fasta", output=OUTPUT)
    output:
        expand("{output}/Binning/{{sample}}/checkm.tsv", output=OUTPUT, sample=set(EXPS['Sample']))
    threads:
        config["threads"]
    params:
        output = lambda wildcards: f'{OUTPUT}/Binning/{wildcards.sample}',
        markerset = config["markerset"],
        iterative = config['do_iterative_binning']
    conda:
        "../envs/binning.yaml"
    script:
        "../scripts/binning.py"

rule dereplication:
    input:
        expand("{output}/Binning/{{sample}}/checkm.tsv", output=OUTPUT)
    output:
        expand("{output}/Binning/{{sample}}/dereplicated/checkm.tsv", output=OUTPUT)
    threads:
        config["threads"]
    params:
        sample = lambda wildcards: wildcards.sample,
        output_dir = lambda wildcards: f'{OUTPUT}/Binning/{wildcards.sample}/dereplicated',
    conda:
        "../envs/binning.yaml"
    shell:
        """
        dRep dereplicate {params.output_dir}/dereplicated -g {params.output_dir}/*.fasta
        checkm lineage_wf -x fasta -r --ali --nt -t {threads} {params.output_dir}/dereplicated --reduced_tree {params.output_dir}/dereplicated --tab_table --file {params.output_dir}/dereplicated/checkm.tsv
        """