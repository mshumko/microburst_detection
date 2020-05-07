import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import numpy as np

import spacepy

from microburst_detection import dirs
from microburst_detection.misc.locate_consecutive_numbers import locateConsecutiveNumbers

class SignalToBackground:
    def __init__(self, counts, cadence, background_width_s):
        """ 
        This class implements the signal to background 
        microburst detection. 

        Warning: Be bareful and feed in the counts array 
        with units of counts/bin, otherwise the Poisson
        significance will be wrong.
        """
        self.counts = counts

        # Check that counts is a DataFrame
        if not isinstance(self.counts, pd.DataFrame):
            self.counts = pd.DataFrame(self.counts)

        self.cadence = cadence
        self.background_width_s = background_width_s
        return

    def significance(self):
        """
        Calculates the number of standard deviations, assuming Poisson
        statistics, that a count value is above a rolling average 
        background of length self.background_width_s.

        Returns a pandas DataFrame object that can be 
        converted to numpy using .to_numpy() method.
        """
        self.rolling_average = self._running_average(self.counts)
        self.n_std = (self.counts-self.rolling_average)/np.sqrt(self.rolling_average+1)
        return self.n_std

    def find_microburst_peaks(self, std_thresh=2):
        """
        This method finds the data intervals where the 
        microburst criteria is satisfied. For for 
        every interval, calculate the time of the highest
        peak.
        """
        self.criteria_idt = np.where(self.n_std >= std_thresh)[0]

        if len(self.criteria_idt) <= 1:
            raise ValueError('No detections found')

        interval_start, interval_end = locateConsecutiveNumbers(self.criteria_idt)
        self.peak_idt = np.nan*np.ones(interval_start.shape[0])

        # Loop over each continous interval and find the peak index.
        for i, (start, end) in enumerate(zip(interval_start, interval_end)):
            if start == end:
                end+=1
            offset = self.criteria_idt[start]
            self.peak_idt[i] = offset + np.argmax(self.counts.loc[self.criteria_idt[start:end]])
        self.peak_idt = self.peak_idt.astype(int)
        return self.peak_idt

    def _running_average(self, counts):
        """
        Calculate the running average of the counts array.
        """
        background_width_samples = int(self.background_width_s/self.cadence)
        return counts.rolling(background_width_samples, center=True).mean()

if __name__ == '__main__':
    hr_name = 'FU4_Hires_2019-09-27_L2.txt'
    background_width_s = 2
    std_thresh = 10
    sc_id = hr_name[2]
    hr_path = pathlib.Path(dirs.firebird_dir(sc_id), hr_name)

    hr = spacepy.datamodel.readJSONheadedASCII(str(hr_path))
    hr['Time'] = pd.to_datetime(hr['Time'])
    cadence = float(hr.attrs['CADENCE'])

    # All of the code to detect microbursts is here.
    s = SignalToBackground(hr['Col_counts'][:, 0], cadence, background_width_s)
    s.significance()
    s.find_microburst_peaks(std_thresh=std_thresh)

    # Now make plots.
    fig, ax = plt.subplots(2, sharex=True)

    ax[0].plot(hr['Time'], hr['Col_counts'][:, 0], 'b', label='Col Counts')
    ax[0].scatter(hr['Time'][s.peak_idt], hr['Col_counts'][s.peak_idt, 0], 
                c='r', marker='*', label='Microburst peaks')
    ax[0].plot(hr['Time'], s.rolling_average, 'r', label=f'{background_width_s} s average')
    ax[1].plot(hr['Time'], s.n_std, label='std above background')
    ax[1].axhline(std_thresh, c='k', label='std thresh')

    ax[0].set(ylabel='FIREBIRD Col counts\n[counts/bin]', 
            title='Signal-to-Background Microburst Detection')
    ax[1].set(ylabel=f'Std above {background_width_s} s\nrunning average',
            xlabel='Time')

    for a in ax:
        a.legend()
    plt.tight_layout()
    plt.show()
