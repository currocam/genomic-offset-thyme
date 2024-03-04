rule misidentified_snps:
    input:
        input=expand("steps/slim/m1.1_s{seed}.vcf", seed=range(500, 550)),
        script="scripts/01_misidentified_snps.jl",
    output:
        expand(
            "steps/misidentified_snps/{basename}.csv",
            basename=[
                "only_causal_snps_one_causal_pred",
                "only_causal_snps_causal_counfounded_pred",
                "all_snps_one_causal_pred",
                "all_snps_two_predictors_confounded",
                "all_snps_two_predictors_uncorrelated",
                "sample_snps_two_predictors_uncorrelated",
                "sample_snps_two_predictors_uncorrelated_unlinked",
            ],
        ),
        touch("steps/misidentified_snps/.done"),
    resources:
        mem_mb=3000,
        runtime=100,
    shell:
        "julia {input.script}"


rule gain_temp:
    input:
        expand("steps/slim/m2.1_s{seed}.vcf", seed=range(500, 551)),
        expand("steps/slim/m2.2_s{seed}.vcf", seed=range(500, 551)),
