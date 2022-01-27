import pathlib
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from microburst_detection import config
from microburst_detection.misc.load_firebird import readJSONheadedASCII
from microburst_detection.misc.locate_consecutive_numbers import locateConsecutiveNumbers

class SignalToBackground:
    def __init__(self, counts, cadence, background_width_s, microburst_width_s):
        """ 
        This class implements the signal to background 
        microburst detection. This method is a generalization 
        of the O'Brien 2003 burst parameter that uses a 0.5 
        second baseline instead of the longer baselines used
        in the examples here.

        To reproduce the O'Brien burst parameter, set 
        background_width_s = 0.5 and microburst_width_s=0.1.

        Parameters
        ----------
        counts : array
            Array of counts. Should be continuous
        cadence : float
            Instrument cadence (seconds)
        microburst_width_s : float
            The duration to seconds integrate the signal, i.e. the n100 
            parameter in the O'Brien 2003 paper.
        background_width_s : float
            The baseline width in seconds to calculate the running mean,
            i.e. the a500 parameter in the O'Brien paper.
        """
        self.counts = counts

        # Check that counts is a DataFrame
        if not isinstance(self.counts, pd.DataFrame):
            self.counts = pd.DataFrame(self.counts)

        self.cadence = cadence
        self.background_width_s = background_width_s
        self.microburst_width_s = microburst_width_s
        return

    def significance(self):
        """
        Calculate the number of background standard deviations, 
        assuming Poisson statistics, that a count value is above
        a rolling average background of length self.background_width_s.

        Returns a pandas DataFrame object that can be 
        converted to numpy using .to_numpy() method.
        """
        self.microburst_rolling_average = self._running_average(self.counts, 
            self.microburst_width_s)
        self.background_rolling_average = self._running_average(self.counts, 
            self.background_width_s)
        
        diff = (self.microburst_rolling_average-self.background_rolling_average)
        self.n_std = diff/np.sqrt(self.background_rolling_average+1)
        return self.n_std

    def find_microburst_peaks(self, std_thresh=2):
        """
        This method finds the data intervals where the 
        microburst criteria is satisfied. Then for
        every interval, calculate the time of the highest
        peak.

        Parameters
        ----------
        std_thresh : float
            The baseline standard deviation threshold above the baseline
            that the data point must be to satisfy the microburst criteria
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

    def _running_average(self, counts, time_window_s):
        """
        Calculate the running average of the counts array.
        """
        background_width_samples = int(time_window_s/self.cadence)
        return counts.rolling(background_width_samples, center=True).mean()


class FirebirdSignalToBackground(SignalToBackground):
    def __init__(self, counts, cadence, background_width_s, microburst_width_s):
        """
        This child class of SignalToBackground finds peaks but reports
        the standard deviations for all 6 FIREBIRD channels.
        """
        super().__init__(counts, cadence, background_width_s, microburst_width_s)
        return

    def significance(self):
        """
        Calculates the number of standard deviations, assuming Poisson
        statistics, that a count value is above a rolling average 
        background of length self.background_width_s. Does this for the
        6 FIREBIRD channels.

        Returns a pandas DataFrame object that can be 
        converted to numpy using .to_numpy() method.
        """
        self.rolling_average = self._running_average(self.counts)
        self.n_std = np.subtract(self.counts, self.rolling_average)/np.sqrt(self.rolling_average+1)
        return self.n_std

    def find_microburst_peaks(self, std_thresh=2, detect_channel=0):
        """
        This method finds the data intervals where the 
        microburst criteria is satisfied. For for 
        every interval, calculate the time of the highest
        peak.
        """
        self.criteria_idt = np.where(self.n_std.loc[:, detect_channel] >= std_thresh)[0]

        if len(self.criteria_idt) <= 1:
            raise ValueError('No detections found')

        interval_start, interval_end = locateConsecutiveNumbers(self.criteria_idt)
        self.peak_idt = np.nan*np.ones(interval_start.shape[0])

        # Loop over each continuos interval and find the peak index.
        for i, (start, end) in enumerate(zip(interval_start, interval_end)):
            if start == end:
                end+=1
            offset = self.criteria_idt[start]
            self.peak_idt[i] = offset + np.argmax(self.counts.loc[self.criteria_idt[start:end], detect_channel])
        self.peak_idt = self.peak_idt.astype(int)
        return self.peak_idt

    def _running_average(self, counts):
        """
        Calculate the running average of the counts array.
        """
        background_width_samples = int(self.background_width_s/self.cadence)
        return counts.rolling(background_width_samples, center=True, axis=0).mean()


if __name__ == '__main__':
    # Detection parameters
    background_width_s = 2
    microburst_width_s = 0.1
    std_thresh = 10

    # Load the HiRes data
    sc_id = 4
    hr_date = datetime(2019, 9, 27) 
    search_str = f'FU{sc_id}*Hires*{hr_date.year}*{hr_date.month}*{hr_date.day}*L2*'
    hr_paths = list(pathlib.Path(config.FB_DIR, ).rglob(search_str))
    assert len(hr_paths) == 1, (f'A unique HiRes path not found.\n'
                                f'hr_paths={hr_paths} search_str={search_str}')
    hr_path = str(hr_paths[0])
    hr = readJSONheadedASCII(str(hr_path))
    hr['Time'] = pd.to_datetime(hr['Time'])
    cadence = float(hr.attrs['CADENCE'])

    # All of the code to detect microbursts is here.
    s = FirebirdSignalToBackground(hr['Col_counts'], cadence, 
                                microburst_width_s, background_width_s)
    s.significance()
    s.find_microburst_peaks(std_thresh=std_thresh)

    # Now make plots.
    fig, ax = plt.subplots(2, sharex=True)
    colors = ['r', 'g', 'b', 'k', 'c', 'y']

    for i in range(6):
        ax[0].plot(hr['Time'], hr['Col_counts'][:, i], colors[i], label=f'Ch {i}')
        ax[0].plot(hr['Time'], s.rolling_average[i], colors[i], ls='--')
        ax[1].plot(hr['Time'], s.n_std[i], colors[i])

    ax[0].scatter(hr['Time'][s.peak_idt], hr['Col_counts'][s.peak_idt, 0], 
                c='r', marker='*', label='Microburst peaks')
    # ax[0].plot(hr['Time'], s.rolling_average, 'r', label=f'{background_width_s} s average')
    
    ax[1].axhline(std_thresh, c='k', label='std thresh')

    ax[0].set(ylabel='FIREBIRD Col counts\n[counts/bin]', 
            title='Signal-to-Background Microburst Detection')
    ax[1].set(ylabel=f'Std above {background_width_s} s\nrunning average',
            xlabel='Time')

    for a in ax:
        a.legend()
    plt.tight_layout()
    plt.show()
