rule plink_to_vcf:
    input:
        plink = [
            f"data/{config['GENOTYPE_BASE_FILE']}.bed",
            f"data/{config['GENOTYPE_BASE_FILE']}.bim",
            f"data/{config['GENOTYPE_BASE_FILE']}.fam"
        ]
    output:
        temp(f"steps/plink_to_vcf/{config['GENOTYPE_BASE_FILE']}.vcf")
    conda:
        "../envs/plink.yaml"
    params:
        plink = f"{config['GENOTYPE_BASE_FILE']}"
    shell:
        """
        plink --bfile data/{params.plink} --recode vcf --out steps/plink_to_vcf/{params.plink} --allow-extra-chr
        """

rule rename_vcf:
    input:
        vcf = f"steps/plink_to_vcf/{config['GENOTYPE_BASE_FILE']}.vcf",
        metadata = config['METADATA_FILE']
    output:
        f"steps/renamed_vcf/{config['GENOTYPE_BASE_FILE']}_renamed.vcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        awk '{{print $1}}' {input.metadata} | tail -n +2 > steps/renamed_vcf/samples_list.txt & bcftools reheader --samples steps/renamed_vcf/samples_list.txt {input.vcf} -o {output}
        """