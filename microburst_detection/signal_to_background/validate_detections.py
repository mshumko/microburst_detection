# Validates the microburst catalogs.
import pathlib
from datetime import date, datetime

import progressbar
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates

from signal_to_background import config
from microburst_detection.misc.load_firebird import readJSONheadedASCII


plot_window_s = 2

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
    return hr

for sc_id in [3,4]:
    for cat_id in [1]:
        catalog_name = f'FU{sc_id}_microburst_catalog_{cat_id:02d}.csv'
        catalog_path = pathlib.Path(config.PROJECT_DIR, 'data', catalog_name)

        save_dir = pathlib.Path(catalog_path.parents[0], 'validation_plots', 
            catalog_path.name.split('.')[0])
        save_dir.mkdir(parents=True, exist_ok=True)

        cat = pd.read_csv(catalog_path, index_col=0, parse_dates=True)

        _, ax = plt.subplots()

        
        current_date = datetime.min
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

            # 20150202_123907_056000_saturated=0_timegap=0_nzeros=0_L=XX_MLT=XX_lat=XX_lon=XX_microburst.png
            save_time = index.strftime("%Y%m%d_%H%M%S_%f")
            save_name = (f'{save_time}_saturated={int(row["saturated"])}_timegap={int(row["time_gap"])}_'
                f'nzeros={int(row["n_zeros"])}_L={round(row["McIlwainL"], 1)}_MLT={round(row["MLT"], 1)}_'
                f'lat={round(row["Lat"], 1)}_lon={round(row["Lon"], 1)}_'
                f'microbursts.png')
            save_path = pathlib.Path(save_dir, save_name)
            plt.savefig(save_path)
            ax.clear()