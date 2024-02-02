from snakemake.shell import shell
import pandas as pd

vcf = snakemake.input.vcf
metadata = snakemake.input.metadata
outfile = snakemake.output[0]

df = pd.read_csv(metadata, sep="\t")
shell(f"mkdir -p {outfile}")
for site in df.site_id.unique():
    # Get samples for this site
    samples = df[df.site_id == site].sample_id.tolist()
    # Filter vcf
    shell(f"bcftools view -s {','.join(samples)} {vcf} > {outfile}/{site}.vcf")
