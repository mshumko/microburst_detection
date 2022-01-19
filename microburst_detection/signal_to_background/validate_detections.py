# Validates the microburst catalogs.
import pathlib
from datetime import date

import progressbar
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates

from signal_to_background import config
from microburst_detection.misc.load_firebird import readJSONheadedASCII

plot_window_s = 10
sc_id = 4
catalog_name = 'FU4_microburst_catalog_02.csv'
catalog_path = pathlib.Path(config.PROJECT_DIR, 'data', catalog_name)

good_save_dir = pathlib.Path(catalog_path.parents[0], 'validation_plots', 
    catalog_path.name.split('.')[0], 'good')
bad_save_dir = pathlib.Path(catalog_path.parents[0], 'validation_plots', 
    catalog_path.name.split('.')[0], 'bad')
good_save_dir.mkdir(parents=True, exist_ok=True)
bad_save_dir.mkdir(parents=True, exist_ok=True)

cat = pd.read_csv(catalog_path, index_col=0, parse_dates=True)

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
    hr = readJSONheadedASCII(hr_paths[0])
    # hr['Time'] = pd.to_datetime(hr['Time'])
    return hr

current_date = date.min

fig, ax = plt.subplots()

for index, row in progressbar.progressbar(cat.iterrows(), max_value=cat.shape[0]):
    if current_date != index.date():
        hr = load_hr(index)
        current_date = index.date()
    dt = pd.Timedelta(seconds=plot_window_s/2)
    time_range = (index-dt, index+dt)

    idt = np.where(
        (hr['Time'] > time_range[0]) &
        (hr['Time'] < time_range[1])
        )[0]
    idt_peak = np.where(hr['Time'] == index)[0]
    ax.plot(hr['Time'][idt], hr['Col_counts'][idt, 0], c='k')
    ax.scatter(hr['Time'][idt_peak], hr['Col_counts'][idt_peak, 0], marker='*', s=200, c='r')

    ax.set(
        xlim=time_range, xlabel='Time', 
        ylabel=f'Counts/{1000*float(hr.attrs["CADENCE"])} ms',
        title=index.strftime("%Y-%m-%d %H:%M:%S.%f\nmicroburst validation")
        )
    s = (
        f'time_gap={row["time_gap"]}\nsaturated={row["saturated"]}\n'
        f'n_zeros={row["n_zeros"]}\n\n'
        f'L={round(row["McIlwainL"], 1)}\n'
        f'MLT={round(row["MLT"], 1)}\n'
        f'(lat,lon)=({round(row["Lat"], 1)}, {round(row["Lon"], 1)})'
        )
    ax.text(0.7, 1, s, va='top', transform=ax.transAxes, color='red')
    locator=matplotlib.ticker.MaxNLocator(nbins=5)
    ax.xaxis.set_major_locator(locator)
    fmt = matplotlib.dates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(fmt)

    plt.tight_layout()
    save_name = index.strftime("%Y%m%d_%H%M%S_%f_microburst.png")
    if (row["time_gap"] == 1) or (row["saturated"] == 1):
        save_path = pathlib.Path(bad_save_dir, save_name)
    else:
        save_path = pathlib.Path(good_save_dir, save_name)
    plt.savefig(save_path)
    ax.clear()