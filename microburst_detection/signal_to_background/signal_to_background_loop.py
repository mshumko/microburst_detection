import numpy as np
import pandas as pd
import pathlib
import progressbar

import spacepy

from microburst_detection.signal_to_background import signal_to_background
from microburst_detection import config

class SignalToBackgroundLoop:
    def __init__(self, sc_id, microburst_width_s, background_width_s, std_thresh, 
        detect_channel=0, catalog_columns=None):
        """
        This program uses signal_to_background detection code to
        loop over all of the FIREBIRD data and detect all 
        microbursts.

        Parameters
        ----------
        sc_id : int
            Spacecraft id. Either 3 or 4
        microburst_width_s: float
            The microburst width to use for the running mean.
        background_width_s : float
            The baseline width in time to calculate the running mean
        std_thresh : float
            The baseline standard deviation threshold above the baseline
            that the data point must be to satisfy the microburst criteria
        detect_channel : int
            The FIREBIRD energy channel number to use for std_thresh 
            critera. This is channel 0 by default.
        catalog_columns : list
            What catalog to save in the catalog. If None, the keys are a 
            combination of HiRes keys, collimated count keys, 
            signal-to-backround keys, and a saturated key.
        """
        self.sc_id = sc_id
        self.microburst_width_s = microburst_width_s
        self.background_width_s = background_width_s
        self.microburst_width_s = microburst_width_s
        self.std_thresh = std_thresh
        self.detect_channel = detect_channel

        if catalog_columns is None:
            self.hr_keys = ['Time', 'Lat', 'Lon', 'Alt', 
                            'McIlwainL', 'MLT', 'kp']
            self.count_keys = [f'counts_s_{i}' for i in range(6)]
            self.sig_keys = [f'sig_{i}' for i in range(6)]
            self.catalog_columns = self.hr_keys + self.count_keys + self.sig_keys + ['saturated']
        else:
            self.catalog_columns = catalog_columns

        # Find all of the HiRes files
        search_str = f'FU{sc_id}_Hires_*L2.txt'
        self.hr_paths = sorted(pathlib.Path(config.FB_DIR).rglob(search_str))
        return

    def loop(self):
        """
        Loop over all the HiRes data and run the signal_to_background
        microburst detector on every day. For the detected microbursts
        save a handful of columns specified by the save_keys kwarg to
        self.microburst_list.
        """        
        self.microburst_list = pd.DataFrame(columns=self.catalog_columns)

        for hr_path in progressbar.progressbar(self.hr_paths, redirect_stdout=True):
            self.hr = spacepy.datamodel.readJSONheadedASCII(str(hr_path))
            self.hr['Time'] = pd.to_datetime(self.hr['Time'])
            self.cadence = float(self.hr.attrs['CADENCE'])
                
            # All of the code to detect microbursts is here.
            self.s = signal_to_background.FirebirdSignalToBackground(
                self.hr['Col_counts'], self.cadence, 
                self.background_width_s, 
                self.microburst_width_s
                )
            self.s.significance()
            try:
                self.s.find_microburst_peaks(std_thresh=self.std_thresh, 
                                        detect_channel=self.detect_channel)
            except ValueError as err:
                if str(err) == 'No detections found':
                    continue
                else:
                    raise

            valid_peaks = self._time_gaps()

            daily_microburst_list = pd.DataFrame(
                data=np.nan*np.ones((len(valid_peaks), len(self.catalog_columns)), dtype=object), 
                columns=self.catalog_columns
                )
            daily_microburst_list.loc[:, self.hr_keys] = np.array(
                [self.hr[col][valid_peaks] for col in self.hr_keys],
                dtype=object
                ).T
            daily_microburst_list.loc[:, self.count_keys] = self.hr['Col_counts'][valid_peaks, :]/self.cadence
            daily_microburst_list.loc[:, self.sig_keys] = self.s.n_std.loc[valid_peaks, :].to_numpy()
                                            
            self.microburst_list = pd.concat((self.microburst_list, daily_microburst_list))
            
        self.microburst_list = self.microburst_list.reset_index()
        del(self.microburst_list['index'])  # Duplicate
        return self.microburst_list

    def save_microbursts(self, save_name=None):
        """
        Save the microburst list to a csv file in the data/ directory.
        If the directory does not exist, one will be created.
        """
        save_dir = pathlib.Path(config.PROJECT_DIR, 'data')

        if not save_dir.is_dir():
            save_dir.mkdir()
            print(f'Made directory at {save_dir}')

        if save_name is None:
            save_name = f'FU{self.sc_id}_microbursts.csv'

        self.microburst_list.to_csv(pathlib.Path(save_dir, save_name), 
                                    index=False)
        return

    def _time_gaps(self, width_s=None, max_time_gap=None):
        """
        For each microburst, check if the time stamps within 
        int(width_s/self.cadence) data points have a time 
        difference less than time_gap. 

        width_s: int
            Used to calculate the time window around each microburst.
            10 seconds if None.
        max_time_gap: float
            The maximum time difference between time stamps within the
            window_s to keep the microburst. 10*self.cadence if None. 
        """
        if width_s is None:
            width_s = 10
        width_dp = width_s//(self.cadence*2)
        if max_time_gap is None:
            max_time_gap = 10*self.cadence

        valid_peaks = np.array([])
        for peak_idt in self.s.peak_idt:
            filtered_times = self.hr['Time'][peak_idt-width_dp:peak_idt+width_dp]
            dt = np.array([(tf-ti).total_seconds() for tf, ti in zip(filtered_times[1:], filtered_times[:-1])])[0]
            if np.max(dt) < max_time_gap:
                valid_peaks = np.append(valid_peaks, peak_idt)
        return valid_peaks
   

if __name__ == '__main__':
    sc_id = 4
    microburst_width_s = 0.1
    background_width_s = 0.5
    std_thresh = 10

    s = SignalToBackgroundLoop(sc_id, microburst_width_s, background_width_s, std_thresh)
    s.loop()
    s.save_microbursts()