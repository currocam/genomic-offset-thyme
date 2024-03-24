# This script assumes that the user has already authenticated with Earth Engine
# You can do this by (a) install gcloud, (b) running `earthengine authenticate` 
# and (c) setting a project via gcloud auth application-default set-quota-project <project_id>

import ee
import pandas as pd
from datetime import datetime

# Trigger the authentication flow.
ee.Authenticate()
# Initialize the library.
ee.Initialize()

# Read metadata
metadata = pd.read_csv('data/genos01_02_2024/thymus_2016_sample_info_clean233samples.tsv', sep='\t')
# Select only site_id  lon_deci  lat_deci
meta_sites = metadata[['site_id', 'lon_deci', 'lat_deci']].drop_duplicates(subset = 'site_id')

# Initial date of interest (inclusive).
i_date = '1978-01-01'
# Final date of interest (exclusive).
f_date = '2025-01-01'
dataset = dataset = ee.ImageCollection("ECMWF/ERA5/MONTHLY").filterDate(i_date, f_date)

def handle_row(site_id, lon_deci, lat_deci):
  point = ee.Geometry.Point(lon_deci, lat_deci)
  raw = dataset.getRegion(point, scale=1).getInfo()
  df = pd.DataFrame(raw[1:len(raw)], columns=raw[0])
  df['site_id'] = site_id
  return df

# Apply the function to each row
res = pd.concat(
  [handle_row(site_id, lon_deci, lat_deci) for site_id, lon_deci, lat_deci in meta_sites.values]
)
res.epoch = res.time.astype(int)
res.time = res.time.apply(lambda epoch: datetime.fromtimestamp(epoch/1000).strftime('%Y-%m-%d'))
res.to_csv('results/ecological/ERA5_233samples.csv', index=False)