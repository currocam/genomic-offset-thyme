rule plink_ld:
    input:
        multiext("steps/slim/{sample}", ".bed", ".bim", ".fam"),
    output:
        "steps/pairwise_ld/{sample}.ld"
    resources:
        mem_mb=300,
        plink_memory=150,
        runtime=5,
    conda:
        "../envs/plink.yaml"        
    log:
        "logs/plink/{sample}_pairwise_ld.log"
    shell:
        """
        plink --memory {resources.plink_memory} \
        --bfile steps/slim/{wildcards.sample} \
        --r2 --matrix \
        --out steps/pairwise_ld/{wildcards.sample} \
        > {log} 2> {log}
        """