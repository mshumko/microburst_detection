import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import numpy as np

import spacepy

from microburst_detection import dirs

class SignalToBackground:
    def __init__(self, cadence, background_width_s):
        """ 
        This class implements the signal to background 
        microburst detection. 

        Warning: Be bareful and feed in the counts array 
        with units of counts/bin, otherwise the Poisson
        significance will be wrong.
        """
        self.cadence = cadence
        self.background_width_s = background_width_s
        return

    def significance(self, counts):
        """
        Calculates the number of standard deviations, assuming Poisson
        statistics, that a count value is above a rolling average 
        background of length self.background_width_s.

        Returns a pandas DataFrame object that can be 
        converted to numpy using .to_numpy() method.
        """
        # Check that counts is a DataFrame
        if not isinstance(counts, pd.DataFrame):
            counts = pd.DataFrame(counts)

        rolling_average = self._running_average(counts)
        self.n_std = (counts-rolling_average)/np.sqrt(rolling_average+1)
        return self.n_std

    def _running_average(self, counts):
        """
        Calculate the running average of the counts array.
        """
        background_width_samples = int(self.background_width_s/self.cadence)
        return counts.rolling(background_width_samples).mean()

if __name__ == '__main__':
    hr_name = 'FU4_Hires_2019-09-27_L2.txt'
    background_width_s = 2
    sc_id = hr_name[2]
    hr_path = pathlib.Path(dirs.firebird_dir(sc_id), hr_name)

    hr = spacepy.datamodel.readJSONheadedASCII(str(hr_path))
    hr['Time'] = pd.to_datetime(hr['Time'])
    cadence = float(hr.attrs['CADENCE'])

    s = SignalToBackground(cadence, background_width_s)
    s.significance(hr['Col_counts'][:, 0])

    fig, ax = plt.subplots(2, sharex=True)

    ax[0].plot(hr['Time'], hr['Col_counts'][:, 0])
    ax[1].plot(hr['Time'], s.n_std)

    plt.show()
