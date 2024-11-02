rule summary_report:
    input:
        expand("{output}/MOSCA_{sample}_General_Report.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx",
        f"{OUTPUT}/DE_analysis/condition_treated_results.tsv" if has_expression_data else [],
        (expand("{output}/Binning/{sample}/checkm.tsv", output=OUTPUT, sample=set(EXPS['Sample']))
         if config['do_binning'] else [])
    output:
        f"{OUTPUT}/MOSCA_Versions_Report.xlsx",
        f"{OUTPUT}/MOSCA_Summary_Report.tsv",
        f"{OUTPUT}/MOSCA_results.zip"
    threads:
        1
    params:
        output=OUTPUT,
        cutoff=config["significance_threshold"],
        has_expression_data=has_expression_data
    conda:
        "../envs/summary.yaml"
    script:
        "../scripts/summary_report.py"
