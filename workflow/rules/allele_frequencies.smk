rule plink_to_vcf:
    input:
        plink = [
            f"data/{config['GENOTYPE_BASE_FILE']}.bed",
            f"data/{config['GENOTYPE_BASE_FILE']}.bim",
            f"data/{config['GENOTYPE_BASE_FILE']}.fam"
        ]
    output:
        f"steps/plink_to_vcf/{config['GENOTYPE_BASE_FILE']}.vcf"
    conda:
        "../envs/plink.yaml"
    params:
        plink = f"{config['GENOTYPE_BASE_FILE']}"
    shell:
        """
        plink --bfile data/{params.plink} --recode vcf --out steps/plink_to_vcf/{params.plink} --allow-extra-chr
        """