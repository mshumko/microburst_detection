import pathlib
import subprocess

import numpy as np
import pandas as pd
import progressbar
import matplotlib.pyplot as plt

from microburst_detection.signal_to_background import signal_to_background
from microburst_detection.misc.load_firebird import readJSONheadedASCII
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
            self.catalog_columns = (self.hr_keys + self.count_keys + 
                self.sig_keys + ['time_gap', 'saturated', 'n_zeros'])
        else:
            self.catalog_columns = catalog_columns

        # Find all of the HiRes files
        search_str = f'FU{sc_id}_Hires_*L2.txt'
        self.hr_paths = sorted(pathlib.Path(config.FB_DIR).rglob(search_str))
        return

    def loop(self, test_plots=False):
        """
        Loop over all the HiRes data and run the signal_to_background
        microburst detector on every day. For the detected microbursts
        save a handful of columns specified by the save_keys kwarg to
        self.microburst_list.
        """        
        self.microburst_list = pd.DataFrame(columns=self.catalog_columns)

        for hr_path in progressbar.progressbar(self.hr_paths, redirect_stdout=True):
            self.hr = readJSONheadedASCII(hr_path)
            self.cadence = self.hr.attrs['CADENCE']
                
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
            
            dropout = self._dropout()
            daily_microburst_list = pd.DataFrame(
                data=np.nan*np.ones((len(self.s.peak_idt), len(self.catalog_columns)), dtype=object), 
                columns=self.catalog_columns
                )
            daily_microburst_list.loc[:, self.hr_keys] = np.array(
                [self.hr[col][self.s.peak_idt] for col in self.hr_keys],
                dtype=object
                ).T
            daily_microburst_list.loc[:, self.count_keys] = self.hr['Col_counts'][self.s.peak_idt, :]/self.cadence
            daily_microburst_list.loc[:, self.sig_keys] = self.s.n_std.loc[self.s.peak_idt, :].to_numpy()
            daily_microburst_list.loc[:, 'time_gap'] = self._time_gaps()
            daily_microburst_list.loc[:, 'saturated'] = dropout[self.s.peak_idt]
            daily_microburst_list.loc[:, 'n_zeros'] = self._number_of_nearby_zeros()
                                            
            self.microburst_list = pd.concat((self.microburst_list, daily_microburst_list))

            if test_plots:
                fig, ax = plt.subplots(2, sharex=True)
                ax[0].plot(self.hr['Time'], self.hr['Col_counts'][:, self.detect_channel], c='k')
                ax[0].scatter(
                    self.hr['Time'][self.s.peak_idt], 
                    self.hr['Col_counts'][self.s.peak_idt, self.detect_channel],
                    marker='X', s=100, c='r', alpha=dropout[self.s.peak_idt]
                    )
                ax[0].scatter(
                    self.hr['Time'][self.s.peak_idt], 
                    self.hr['Col_counts'][self.s.peak_idt, self.detect_channel],
                    marker='*', s=200, c='r', alpha=1-dropout[self.s.peak_idt]
                    )
                ax[1].plot(self.hr['Time'], dropout)
                # ax[1].plot(self.hr['Time'], dropout)
                plt.show()
            
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
            counter = 0
            while True:
                save_name = f'FU{self.sc_id}_microburst_catalog_{counter:02d}.csv'.format(counter)
                save_path = pathlib.Path(save_dir, save_name)
                if not save_path.exists():
                    break
                counter += 1
        else:
            save_path = pathlib.Path(save_dir, save_name)

        self.microburst_list.to_csv(save_path, index=False)
        self._save_log(save_path)
        return

    def _save_log(self, save_path):
        """
        Save the catalog log
        """
        log_path = pathlib.Path(save_path.parents[0], 'catalog_log.csv')

        # Log the saved catalog info.
        git_revision_hash = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD']
            ).strip().decode()
        log = pd.DataFrame(
            index=[0],
            data={ 
                'time':pd.Timestamp.today(),
                'catalog_name':save_path.name,
                'burst_params':repr(self),
                'git_revision_hash':git_revision_hash
                })
        # Determine if the header needs to be written
        if log_path.exists():
            header=False
        else:
            header=True
        log.to_csv(log_path, 
                mode='a', header=header, index=False)
        return

    def _time_gaps(self, width_s=5, max_time_gap=None):
        """
        For each microburst, check if the time stamps within 
        int(width_s/self.cadence) data points have a time 
        difference less than time_gap. 

        width_s: int
            Used to identify the time window (in data points) 
            around each microburst to test for time gaps.
        max_time_gap: float
            The maximum allowable time difference between time 
            stamps within the window_s okay keep the microburst. 
            10*self.cadence if None. 
        """
        width_dp = int(width_s/(self.cadence*2))
        if max_time_gap is None:
            max_time_gap = 5*self.cadence

        near_gap = np.ones_like(self.s.peak_idt)  # Default to all near a time gap.
        for i, peak_idt in enumerate(self.s.peak_idt):
            # A time gap if the detection was made in the very begining or end of the day.
            if (peak_idt-width_dp < 0) or (peak_idt+width_dp >= len(self.hr['Time'])):
                continue
            filtered_times = self.hr['Time'][peak_idt-width_dp:peak_idt+width_dp]
            dt = np.array([(tf-ti).total_seconds() for tf, ti in zip(filtered_times[1:], filtered_times[:-1])])
            if np.max(dt) < max_time_gap:
                near_gap[i] = 0
        return near_gap

    def _number_of_nearby_zeros(self, width_s=5):
        """
        Look for, and count how many zero counts were near each microburst.

        width_s: int
            Used to identify the time window (in data points) 
            around each microburst to look for zeros.
        """
        width_dp = int(width_s/(self.cadence*2))
        
        n_zeros = np.zeros_like(self.s.peak_idt)  # Default to all near a time gap.
        for i, peak_idt in enumerate(self.s.peak_idt):
            if (peak_idt-width_dp < 0):
                start_index = 0 
            else: 
                start_index = peak_idt-width_dp

            if (peak_idt+width_dp >= len(self.hr['Time'])):
                end_index = len(self.hr['Time'])
            else:
                end_index = peak_idt+width_dp
                
            n_zeros[i] = sum(self.hr['Col_counts'][start_index:end_index, 0]==0)
        return n_zeros


    def _dropout(self, derivative_thresh=300, quarantine_dp=20):
        """
        Identify dropouts in the FIREBIRD data using the derivative approach.

        derivative_thresh: float
            How much the counts need to change by (increase or decrease) over one data
            point.
        quarantine_dp: int
            How many data points around the dropout to flag as affected by the dropout.
        """
        counts = self.hr['Col_counts'][:, self.detect_channel]
        # Techincally there is a division here, but we can ignore it since were working in count space.
        dc_dt = counts[1:] - counts[:-1]
        dropouts = np.zeros_like(counts, dtype=int)

        for i, (dc_dt_i, dc_dt_f) in enumerate(zip(dc_dt[:-1], dc_dt[1:])):
            # Check for a large drop, immediately followed by a large increase 
            # by greather than derivative_thresh counts (per bin). 
            if (dc_dt_i < -derivative_thresh) and (dc_dt_f > derivative_thresh):
                start_index = max(0, i-quarantine_dp)
                end_index = min(i+quarantine_dp, len(dropouts)-1)
                # Mark all data points within +/- quarantine_dp as a dropout.
                dropouts[start_index:end_index] = 1
        return dropouts

    def __repr__(self):
        params = (
                f'sc_id={self.sc_id}, '
                f'microburst_width_s={self.microburst_width_s}, '
                f'background_width_s={self.background_width_s}, '
                f'std_thresh={self.std_thresh}'
                )
        return f'{self.__class__.__qualname__}(' + params + ')'
   

if __name__ == '__main__':
    sc_id = 4
    microburst_width_s = 0.1
    background_width_s = 0.5
    std_thresh = 10

    s = SignalToBackgroundLoop(sc_id, microburst_width_s, background_width_s, std_thresh)
    s.loop()
    s.save_microbursts()