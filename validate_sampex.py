### This code loads in the SAMPEX HILT data, and runs Paul's and 
### wavelet detection methods

import numpy as np
import csv
import os
#import sys
from datetime import datetime, timedelta

import burst_parameter.burst_parameter as burst_parameter
import wavelets.waveletAnalysis as waveletAnalysis

import matplotlib.pyplot as plt

class ValidateBurstAlgorithm():
    def __init__(self, fPath):
        self._load_sampex_hilt(fPath)
        return

    def _load_sampex_hilt(self, path):
        """
        State 4: Second 20-msec SSD configuration: 1996-220 thru 2004-182. 
        (Sum is total of SSD1 through SSD4 rates).
        Rate1 : Sum from Time to Time + 20 msec
        Rate2 : Sum from Time + 20 msec to Time + 40 msec
        Rate3 : Sum from Time + 40 msec to Time + 60 msec
        Rate4 : Sum from Time + 60 msec to Time + 80 msec
        Rate5 : SSD4 from Time to Time + 100 msec
        Rate6 : Sum from Time + 80 msec to Time + 100 msec
        """
        with open(path) as f:
            r = csv.reader(f, delimiter=' ')
            keys = next(r)
            rawData = np.array(list(r)).astype(float)
        # Find the start date from filename.
        rdate = path.split('.')[0][-7:]
        startTime = datetime.strptime(rdate, '%Y%j')
        # Populate data dict
        self.helt = {}
        for (i, key) in enumerate(keys):
            if 'time' in key.lower():
                self.helt['dateTime'] = np.array([startTime+timedelta(seconds=dt) 
                                        for dt in rawData[:, i]])
            else:
                self.helt[key] = rawData[:, i]/100#*10/15 # Convert to flux
        return

    def run_burst_param(self, ch='Rate5', thresh=10, N_WIDTH=0.1, A_WIDTH=0.5):
        """

        """
        self.flag = burst_parameter.obrien_burst_param(self.helt[ch], 0.1,
                 N_WIDTH=N_WIDTH, A_WIDTH=A_WIDTH)
        self.detection = np.where(self.flag > thresh)[0]
        return

    def run_wavelet_det(self, thresh=0.1):
        self.w = waveletAnalysis.WaveletDetector(self.helt['Rate5'],
                    self.helt['dateTime'], 0.1)
        self.w.waveletTransform()    
        self.w.waveletFilter(self.w.s0, 1)
        self.w.degenerateInvWaveletTransform()
        self.wDetection = np.where(self.w.dataFlt > thresh)[0]
        return

if __name__ == '__main__':
    fdir = '/home/mike/research/sampex/hilt'
    fname = 'hhrr1999229.txt'
    v = ValidateBurstAlgorithm(os.path.join(fdir, fname))
    v.run_burst_param()
    v.run_wavelet_det()

    ### VISUALIZE DETECTIONS ###
    plt.plot(v.helt['dateTime'], v.helt['Rate5'], 'k')
    plt.scatter(v.helt['dateTime'][v.detection], 
                v.helt['Rate5'][v.detection], 
                marker='x', c='r', s=100)
    plt.scatter(v.helt['dateTime'][v.wDetection], 
                v.helt['Rate5'][v.wDetection], 
                marker='x', c='g', s=100)

    plt.xlim(datetime(1999, 8, 17, 4, 11), datetime(1999, 8, 17, 4, 16))
    #plt.ylim(1, 10**4.5)
    plt.yscale('log')
    plt.xlabel('UTC')
    #plt.ylabel('SAMPEX SSD4 flux (rate 5)')
    #plt.ylabel('SAMPEX SSD4 counts/100 ms/100')
    #plt.ylabel('SAMPEX SSD4 artificial')
    plt.title("O'brien burst parameter validation")
    plt.show()