# Validates the microburst catalogs.
import pathlib
from datetime import date


import progressbar
import pandas as pd
import numpy as np
import spacepy
import matplotlib.pyplot as plt

from signal_to_background import config

plot_window_s = 10
sc_id = 4
catalog_name = 'FU4_microburst_catalog_00.csv'
catalog_path = pathlib.Path(config.PROJECT_DIR, 'data', catalog_name)

save_dir = pathlib.Path(catalog_path.parents[0], 'validation_plots')
save_dir.mkdir(parents=True, exist_ok=True)

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

current_date = date.min

for index, row in progressbar.progressbar(cat.iterrows(), max_value=cat.shape[0]):
    if current_date != index.date():
        hr = load_hr(index)
        current_date = index.date()
    dt = pd.Timedelta(seconds=plot_window_s/2)
    time_range = (index-dt, index+dt)

    _, ax = plt.subplots()

    idt = np.where(
        (hr['Time'] > time_range[0]) &
        (hr['Time'] < time_range[1])
        )[0]
    idt_peak = np.where(hr['Time'] == index)[0]
    ax.plot(hr['Time'][idt], hr['Col_counts'][idt, 0], c='k')
    ax.scatter(hr['Time'][idt_peak], hr['Col_counts'][idt_peak, 0], marker='*', s=200, c='r')

    ax.set(
        xlim=time_range, xlabel='Time', 
        ylabel=f'Counts/{float(hr.attrs["CADENCE"])}',
        title=index.strftime("%Y-%m-%d %H:%M:%S microburst validation")
        )
    ax.text(0.8, 1, f'time_gap={row["time_gap"]}\nsaturated={row["saturated"]}', 
        va='top', transform=ax.transAxes)

    plt.tight_layout()
    save_name = index.strftime("%Y%m%d_%H%M%S_microburst.png")
    save_path = pathlib.Path(save_dir, save_name)
    plt.savefig(save_path)
    del(ax)