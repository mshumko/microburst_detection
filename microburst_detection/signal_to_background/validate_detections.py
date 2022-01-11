# Validates the microburst catalogs.
import pathlib

import pandas as pd
import matplotlib.pyplot as plt

from signal_to_background import config

catalog_name = 'FU4_microburst_catalog_00.csv'
catalog_path = pathlib.Path(config.PROJECT_DIR, 'data', catalog_name)

cat = pd.read_csv(catalog_path, index_col=0, parse_dates=True)
cat = cat[(cat['time_gap'] == 1) | (cat['saturated'] == 1)]
print(cat.shape)
pass