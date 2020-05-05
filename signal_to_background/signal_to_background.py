import pandas as pd

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
        """
        rolling_average = self._running_average(counts)
        self.n_std = (counts-rolling_average)/np.sqrt(rolling_average+1)
        return self.n_std

    def _running_average(self, counts):
        """
        Calculate the running average of the counts array.

        Returns a pandas DataFrame object that can be 
        converted to numpy using .to_numpy() method.
        """
        # Check that counts is a DataFrame
        if not isinstance(counts, pd.DataFrame):
            counts = pd.DataFrame(counts)

        background_width_samples = int(self.background_width_s/self.cadence)
        return counts.rolling(background_width_samples)
