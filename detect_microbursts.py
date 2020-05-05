import numpy as np

import wavelets.waveletAnalysis as waveletAnalysis
import burst_parameter.burst_parameter as burst_parameter
import signal_to_background.signal_to_background as signal_to_background

class Detect:
    def __init__(self, method, **kwargs):

        return

    def detect(self, std_thresh:float=2, corr_thresh:float=None) -> typing.List[int]:
        """
        After running the baseline_significance() and/or rolling_correlation() methods,
        use this method to find where both conditions are true (intersection).
        """
        self.std_thresh = std_thresh
        self.corr_thresh = corr_thresh

        # Find indicies in the AC6A and B data that are significant 
        # above the background
        idx_signif = np.where((self.n_std_a > std_thresh) & 
                            (self.n_std_b > std_thresh))[0]
        if corr_thresh is not None:
            # Find indicies where the temporal time series was 
            # highly correlated
            idx_corr = np.where(self.corr > corr_thresh)[0]
            idx_signif = list(set(idx_signif).intersection(idx_corr))
        # Now find good quality data.
        valid_data = self.valid_data_flag()
        intersect = set(idx_signif).intersection(valid_data)
        self.detections = np.array(list(intersect))
        self.detections = np.sort(self.detections)
        if len(self.detections) <= 1:
            raise ValueError('No detections were found on this day.')
        return self.detections

    def find_peaks(self):
        """
        Given the self.detections array, for each continious interval of 
        indicies find the index with highest count rates in that ineterval 
        from both spacecraft.
        """
        startInd, endInd = locate_consecutive_numbers.locateConsecutiveNumbers(
            self.detections) # Find consecutive numbers to get a max of first
        self.peaks_A = np.nan*np.ones(len(startInd), dtype=int)
        self.peaks_B = np.nan*np.ones(len(startInd), dtype=int)
        # Loop over every microburst detection region (consecutive microburst indicies)
        for i, (st, et) in enumerate(zip(startInd, endInd)):
            if st == et: 
                # If the interval is just one point 
                et += 1
            # Find the max and reindex.
            offset = self.df_a.index[self.detections[st]]
            self.peaks_A[i] = np.argmax(
                    self.df_a.loc[self.detections[st:et], 'dos1rate']) + offset
            self.peaks_B[i] = np.argmax(
                    self.df_b.loc[self.detections[st:et], 'dos1rate']) + offset
        if len(self.peaks_A) <= 1:
            raise ValueError('No detections were found on this day.')
        return

