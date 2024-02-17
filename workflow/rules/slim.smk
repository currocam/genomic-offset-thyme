rule vcf_slim:
    input:
        "src/slim/{model}.slim"
    output:
        "steps/slim/{model}_s{seed}.vcf",
        "steps/slim/{model}_s{seed}.txt"
    resources:
        mem_mb=300,
        runtime=100
    conda:
        "../envs/slim.yaml"
    log:
        "logs/slim/{model}_s{seed}.log"
    shell:
        """
        slim -s {wildcards.seed} -d "outvcf='{output[0]}'" -d "outfile='{output[1]}'" < {input} > {log}
        """
