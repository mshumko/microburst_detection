"""
Wavelet test script to make detections with the FIREBIRD data.
"""
import pathlib
import string
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pylab as plt

from microburst_detection.wavelets import wavelet_analysis
from microburst_detection import config
from microburst_detection.misc.load_firebird import readJSONheadedASCII

plt.rcParams.update({'font.size': 13})

# Load example HiRes data
sc_id = 4
hr_date = datetime(2019, 9, 27) 
search_str = f'FU{sc_id}*Hires*{hr_date.year}*{hr_date.month}*{hr_date.day}*L2*'
hr_paths = list(pathlib.Path(config.FB_DIR, ).rglob(search_str))
assert len(hr_paths) == 1, (f'A unique HiRes path not found.\n'
                            f'hr_paths={hr_paths} search_str={search_str}')
hr_path = str(hr_paths[0])
hr = readJSONheadedASCII(hr_path)

# Convert time array
hr['Time'] = pd.to_datetime(hr['Time'])
cadence = float(hr.attrs['CADENCE'])
cadence_int = int(cadence*1000)

# Filter the data by time to reduce computational time.
if True:
    startTime = datetime(2019, 9, 27, 19, 30, 30)
    endTime = datetime(2019, 9, 27, 19, 31, 5)
    validInd = np.where((hr['Time'] > startTime) & (hr['Time'] < endTime))[0]
    assert len(validInd) > 0, "No data found in the time range specified!"
else:
    validInd = np.arange(hr['Time'].shape[0])

channel=0
count_thresh=0.1
max_width = 1

# Run the wavelet code
waveDet = wavelet_analysis.WaveletDetector(
    hr['Col_counts'][validInd, channel],
    hr['Time'][validInd], 
    cadence=cadence, 
    siglvl=0.95, 
    run_scipt=False, 
    j1=40
    )

waveDet.waveletTransform()
waveDet.waveletFilter(waveDet.s0, max_width)
waveDet.degenerateInvWaveletTransform()
waveDet.TestForMicrobursts(COUNT_THRESH=count_thresh)
waveDet.findMicroburstPeaks()

fig, ax = plt.subplots(3, figsize=(8, 8), sharex=True)

## Plot raw count data, peaks, and wavelet-based significant (microburst) counts
ax[0].plot(hr['Time'][validInd], 
            hr['Col_counts'][validInd, channel],
             'k')
# ax[0].plot(hr['Time'][validInd[waveDet.indicies]],
#             hr['Col_counts'][validInd[waveDet.indicies], channel], 
#             'go', label='Microburst condition satisfied')
ax[0].plot(hr['Time'][validInd[waveDet.peaks]],
            hr['Col_counts'][validInd[waveDet.peaks], channel], 
            'r*', ls='None', markersize=10, label='microburst peak')

ax[2].plot(hr['Time'][validInd], 
            waveDet.dataFlt,
             'k')
ax[2].axhline(count_thresh, ls='--')

ax[0].set(ylabel=f'original counts\n[counts/{cadence_int} ms]', 
        title=f'FIREBIRD unit {sc_id} | collimated detector | {hr_date.date()}')
ax[1].set(ylabel='period [s]')
ax[2].set(ylabel=f'filtered counts', xlabel='Time')
# ax[0].legend()

# Plot wavelet power
waveDet.plotPower(ax[1])
ax[1].axhline(max_width, c='w', lw=1, ls='--')
ax[1].set_ylim(4, None)
# Hatch the periods outside the interval we care about
ax[1].fill_between(hr['Time'][validInd], max_width*np.ones_like(validInd), 10, facecolor="none", hatch="X", edgecolor="w", linewidth=0.0)

subplot_titles = ['original data', 'wavelet power spectrum', 'filtered data']
subplot_text_y_pos = [1, 0.1, 1]
subplot_text_color=['k', 'k', 'k']
for i, (ax_i, title_i, y_pos, c_i) in enumerate(zip(ax, subplot_titles, subplot_text_y_pos, subplot_text_color)):
    ax_i.text(0, y_pos, f'({string.ascii_lowercase[i]}) {title_i}', va='top',
            transform=ax_i.transAxes, fontsize=15, weight='bold', c=c_i)

plt.tight_layout()
plt.show()