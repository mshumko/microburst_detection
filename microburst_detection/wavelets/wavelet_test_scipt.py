"""
Wavelet test script to make detections with the FIREBIRD data.
"""
import pathlib
import string
from datetime import datetime

import spacepy.datamodel
import pandas as pd
import numpy as np
import matplotlib.pylab as plt

from microburst_detection.wavelets import wavelet_analysis
from microburst_detection import config

# Load example HiRes data
sc_id = 4
hr_date = datetime(2019, 9, 27) 
search_str = f'FU{sc_id}*Hires*{hr_date.year}*{hr_date.month}*{hr_date.day}*L2*'
hr_paths = list(pathlib.Path(config.FB_DIR, ).rglob(search_str))
assert len(hr_paths) == 1, (f'A unique HiRes path not found.\n'
                            f'hr_paths={hr_paths} search_str={search_str}')
hr_path = str(hr_paths[0])
hr = spacepy.datamodel.readJSONheadedASCII(hr_path)

# Convert time array
hr['Time'] = pd.to_datetime(hr['Time'])

# Filter the data by time to reduce computational time.
if True:
    startTime = datetime(2019, 9, 27, 18, 26, 30)
    endTime = datetime(2019, 9, 27, 18, 29, 30)
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
    cadence=float(hr.attrs['CADENCE']), 
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

ax[1].plot(hr['Time'][validInd], 
            waveDet.dataFlt,
             'k')
ax[1].axhline(count_thresh, ls='--')

ax[0].set(ylabel='counts [counts/bin]', 
        title=f'FIREBIRD unit {sc_id} | collimated detector | {hr_date.date()}')
ax[1].set(ylabel='counts [counts/bin]')
ax[2].set(xlabel='Time', ylabel='period [s]')
# ax[0].legend()

# Plot wavelet power
waveDet.plotPower(ax[-1])
ax[-1].axhline(max_width, c='w', lw=3)

subplot_titles = ['original data', 'filtered data', 'wavelet power spectrum']
for i, (ax_i, title_i) in enumerate(zip(ax, subplot_titles)):
    ax_i.text(0, 1, f'({string.ascii_lowercase[i]}) {title_i}', va='top',
            transform=ax_i.transAxes, fontsize=15, weight='bold')

plt.tight_layout()
plt.show()