# wavelet test script
import spacepy.datamodel
import pandas as pd
# import spacepy.time
import datetime
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import pathlib

import waveletAnalysis
import dirs

# Load example HiRes data
hr_name = 'FU4_Hires_2019-09-27_L2.txt'
sc_id = 4    
hr_path = str(pathlib.Path(dirs.firebird_dir(sc_id), hr_name))
hr = spacepy.datamodel.readJSONheadedASCII(hr_path)
# Convert time array
hr['Time'] = pd.to_datetime(hr['Time'])

# # Filter the data by time to reduce computational time.
if True:
    startTime = datetime.datetime(2019, 9, 27, 18, 26, 30)
    endTime = datetime.datetime(2019, 9, 27, 18, 29, 30)
    validInd = np.where((hr['Time'] > startTime) & (hr['Time'] < endTime))[0]
    assert len(validInd) > 0, "No data found in the time range specified!"
else:
    validInd = np.arange(hr['Time'].shape[0])

channel = 0

# Run the wavelet code
waveDet = waveletAnalysis.WaveletDetector(
    hr['Col_counts'][validInd, channel],
    hr['Time'][validInd], 
    cadence=float(hr.attrs['CADENCE']), 
    siglvl=0.95, 
    run_scipt=False, 
    j1=40
    )

waveDet.waveletTransform()
waveDet.waveletFilter(waveDet.s0, 10.0)
waveDet.degenerateInvWaveletTransform()
waveDet.TestForMicrobursts(COUNT_THRESH=0.01)
waveDet.findMicroburstPeaks()

fig, ax = plt.subplots(3, figsize=(8, 8), sharex=True)

## Plot raw count data, peaks, and wavelet-based significant (microburst) counts
ax[0].plot(hr['Time'][validInd], 
            hr['Col_counts'][validInd, channel],
             'b')
ax[1].plot(hr['Time'][validInd], 
            waveDet.dataFlt,
             'b')
ax[0].plot(hr['Time'][validInd[waveDet.indicies]],
            hr['Col_counts'][validInd[waveDet.indicies], channel], 
            'go', label='Microburst condition satisfied')
ax[0].plot(hr['Time'][validInd[waveDet.peaks]],
            hr['Col_counts'][validInd[waveDet.peaks], channel], 
            'r*', ls='None', markersize=10, label='microburst peak')

ax[0].set(ylabel='Col counts [counts/bin]', title='Original data')
ax[1].set(ylabel='Col counts [counts/bin]', title='Filtered data')
ax[2].set(xlabel='Time', ylabel='period [s]', title='Wavelet power spectrum')
ax[0].legend()

# Plot wavelet power
waveDet.plotPower(ax[-1])
plt.tight_layout()
plt.show()