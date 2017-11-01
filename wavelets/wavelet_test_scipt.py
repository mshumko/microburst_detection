# wavelet test script

import spacepy.datamodel
import spacepy.time
import datetime
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0, '/home/mike/FIREBIRD/microburst_detector/')
import flag_dropouts
from waveletAnalysis import WaveletDetector

if 'hr' not in globals():
    fname = 'FU4_Hires_2015-02-02_L2.txt'
    fdir = '/home/mike/FIREBIRD/Datafiles/FU_4/hires/level2/'
    
    hr = spacepy.datamodel.readJSONheadedASCII(fdir + fname)
    
    times = spacepy.time.Ticktock(hr['Time']).UTC
    
#startTime = datetime.datetime(2015, 2, 2, 5, 0, 0)
#endTime = datetime.datetime(2015, 2, 2, 7, 0, 0)
startTime = datetime.datetime(2015, 2, 2, 0, 0, 0)
endTime = datetime.datetime(2015, 2, 3, 0, 0, 0)
    
channel = 0
validInd = np.where((times > startTime) & (times < endTime))[0]

assert len(validInd) > 0, "No data found in the time range specified!"

# Run the wavelet code
waveDet = WaveletDetector(hr['Col_counts'][validInd, channel], \
times[validInd], cadence = 18.75E-3, siglvl = 0.7, run_scipt = False)
waveDet.waveletTransform()
waveDet.waveletFilter(waveDet.s0, 1.0)
waveDet.degenerateInvWaveletTransform()
waveDet.TestForMicrobursts(COUNT_THRESH = 0.0)
waveDet.findMicroburstPeaks()
waveDet.savePeakIndicies('/home/mike/FIREBIRD/microburst_detector/wavelets/'+\
'microburst_peak_output/', fname[0:3] + '_' + fname[10:20] + '_peaks')

## Calculate the dropout flag.
#dropoutFlag = flag_dropouts.dropOutFlag(hr['Col_counts'][validInd, channel])
#
## Calculate the smoothed data.
#kernelN = 5
#dataSmoothed = np.convolve(np.ones(kernelN)/kernelN, \
#hr['Col_counts'][validInd, channel], mode = 'same')
#
## Set up plots
#fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'grey')
#gs = gridspec.GridSpec(2, 1)
#dataPlt = fig.add_subplot(gs[0, 0])
#wavePlt = fig.add_subplot(gs[1, 0], sharex = dataPlt)
#flagPlt = dataPlt.twinx()
#
## Plot raw count data, peaks, and wavelet-based significant (microburst) counts
#dataPlt.plot(times[validInd], hr['Col_counts'][validInd, channel], 'b')
#dataPlt.plot(times[validInd], dataSmoothed, 'k')
#dataPlt.plot(times[validInd[waveDet.indicies]], \
#hr['Col_counts'][validInd[waveDet.indicies], channel], 'go')
#dataPlt.plot(times[validInd[waveDet.peaks]], \
#hr['Col_counts'][validInd[waveDet.peaks], channel], 'g*', ls = 'None', \
#markersize = 20)
## Plot the bad data flag/
#flagPlt.plot(times[validInd], dropoutFlag, 'r')
#flagPlt.set_ylim([0, 1.5])
## Plot wavelet power
#waveDet.plotPower(wavePlt)
#gs.tight_layout(fig)