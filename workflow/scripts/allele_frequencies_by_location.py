import os
import pandas as pd
from pysam import VariantFile

indir = snakemake.input[0]
outfile = snakemake.output[0]

# For every .vcf file in the input directory
vcfs = [f for f in os.listdir(indir) if f.endswith(".vcf")]
sites = [f.replace(".vcf", "") for f in vcfs]

df = None
for file, site in zip(vcfs, sites):
    if df is None:
        ids = [rec.id for rec in VariantFile(f"{indir}/{file}") if rec.id]
        maf = [
            rec.info["AC"][0] / rec.info["AN"] for rec in VariantFile(f"{indir}/{file}")
        ]
        df = pd.DataFrame({"id": ids, site: maf})
    else:
        assert df.id.tolist() == [rec.id for rec in VariantFile(f"{indir}/{file}")]
        df[site] = [
            rec.info["AC"][0] / rec.info["AN"] for rec in VariantFile(f"{indir}/{file}")
        ]
df.to_csv(outfile, sep="\t", index=False)
