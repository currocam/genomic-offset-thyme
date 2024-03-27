rule vcf_slim:
    input:
        "src/slim/{model}.slim",
    output:
        "steps/slim/{model}_s{seed}.vcf",
        "steps/slim/{model}_s{seed}.txt",
    wildcard_constraints:
        seed="\d+",
    resources:
        mem_mb=300,
        runtime=100,
    conda:
        "../envs/slim.yaml"
    log:
        "logs/slim/{model}_s{seed}.log",
    shell:
        """
        slim -s {wildcards.seed} -d "outvcf='{output[0]}'" -d "outfile='{output[1]}'" < {input} > {log}
        """

rule slim1:
    input:
        "src/slim/{model}.slim",
    output:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}.vcf",
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}.txt",
    resources:
        mem_mb=300,
        runtime=100,
    wildcard_constraints:
        seed="\d+",
    conda:
        "../envs/slim.yaml"
    log:
        "logs/slim/{model}_s{seed}_nQTL1s{nQTL1s}.log",
    shell:
        """
        slim -s {wildcards.seed} -d "nQTLs1='{wildcards.nQTL1s}'" -d "outvcf='{output[0]}'" -d "outfile='{output[1]}'" < {input} > {log}
        """


rule slim2:
    input:
        "src/slim/{model}.slim",
    output:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.vcf",
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.txt",
    resources:
        mem_mb=300,
        runtime=100,
    wildcard_constraints:
        seed="\d+",
    conda:
        "../envs/slim.yaml"
    log:
        "logs/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.log",
    shell:
        """
        slim -s {wildcards.seed} -d "nQTLs1='{wildcards.nQTL1s}'" -d "nQTLs2='{wildcards.nQTL2s}'" -d "outvcf='{output[0]}'" -d "outfile='{output[1]}'" < {input} > {log}
        """

rule parse_slim1:
    input:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}.vcf",
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}.txt",
    output:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}.Rds",
    resources:
        mem_mb=800,
        runtime=5,
    params:
        script="workflow/scripts/parse_slim.jl",
    shell:
        """
        julia  --project=. {params.script} {input} {output}
        """

rule parse_slim2:
    input:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.vcf",
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.txt",
    output:
        "steps/slim/{model}_s{seed}_nQTL1s{nQTL1s}_nQTL2s_{nQTL2s}.Rds",
    resources:
        mem_mb=800,
        runtime=5,
    params:
        script="workflow/scripts/parse_slim.jl",
    shell:
        """
        julia  --project=. {params.script} {input} {output}
        """

rule vcf2plink:
    input:
        "steps/slim/{sample}.vcf",
    output:
        temp(multiext("steps/slim/{sample}", ".bed", ".bim", ".fam")),
    resources:
        mem_mb=300,
        plink_memory=150,
        runtime=5,
    conda:
        "../envs/plink.yaml"        
    log:
        "logs/plink/{sample}_vct2plink.log"
    shell:
        """
        plink --memory {resources.plink_memory} --vcf {input} --out steps/slim/{wildcards.sample} > {log}
        """
