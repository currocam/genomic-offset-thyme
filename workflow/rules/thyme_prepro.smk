rule go_preprocess:
    input:
       prunned = "data/genos27_03_2024/tVulgragtag_thyme2016CoreSampling_maf0.1_geno0.05_imiss0.1_r20.1_prunned_hetremove.bed",
       unprunned = "data/genos27_03_2024/tVulgragtag_thyme2016CoreSampling_maf0.1_geno0.05_imiss0.1_hetremove.bed",
       annotation_tsv = "data/genos27_03_2024/tVulgragtag_thyme2016CoreSampling_maf0.1_geno0.05_imiss0.1_hetremove.snps_annotation.tsv",
       sample_tsv = "data/genos27_03_2024/thymus_2016_sample_info_clean231samples_hetremove.tsv",
       environment_rds = "data/env_00_12_2019/EnvVariables_5PCs_dec2019.rds"
    output:
        "steps/genomic_offset/genos27_03_2024_prepro.Rds",
    log:
        "logs/genomic_offset/genos27_03_2024_prepro.log",
    resources:
        mem_mb=3000,
        runtime=10,
    script:
        "../scripts/prepro_thyme.jl"


rule ee_ERA5:
    input:
        metadata = "data/genos27_03_2024/thymus_2016_sample_info_clean231samples_hetremove.tsv",
        script = "scripts/05_earth_engine_ERA5.py"
    output:
        "results/ecological/ERA5_231samples.csv",
    log:
        "logs/ecological/ERA5_231samples.log",
    resources:
        mem_mb=1000,
        runtime=5,
    shell:
        "python {input.script} > {log} 2> {log}"

