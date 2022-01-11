# Validates the microburst catalogs.
import pathlib
from datetime import datetime

import progressbar
import pandas as pd
import spacepy
import matplotlib.pyplot as plt

from signal_to_background import config

sc_id = 4
catalog_name = 'FU4_microburst_catalog_00.csv'
catalog_path = pathlib.Path(config.PROJECT_DIR, 'data', catalog_name)

cat = pd.read_csv(catalog_path, index_col=0, parse_dates=True)
cat = cat[(cat['time_gap'] == 1) | (cat['saturated'] == 1)]

def load_hr(date):
    """
    Load the HiRes data.
    """
    search_str = f'FU{sc_id}_Hires_{date.strftime("%Y-%m-%d")}_L2.txt'
    hr_paths = sorted(pathlib.Path(config.FB_DIR).rglob(search_str))
    assert len(hr_paths) == 1, (
        f'{len(hr_paths)} HiRes paths found in {pathlib.Path(config.FB_DIR)} '
        f'that match {search_str}.'
        )
    hr = spacepy.datamodel.readJSONheadedASCII(str(hr_paths[0]))
    hr['Time'] = pd.to_datetime(hr['Time'])
    return hr

current_date = datetime.min

for index, row in progressbar.progressbar(cat.iterrows(), max_value=cat.shape[0]):
    if current_date.date() != index.date():
        hr = load_hr(index)
        current_date = index.date()
    pass