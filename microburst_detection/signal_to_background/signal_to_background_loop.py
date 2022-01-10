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
        self.microburst_list = pd.DataFrame(columns=self.save_keys)

        for hr_path in progressbar.progressbar(self.hr_paths, redirect_stdout=True):
            hr = spacepy.datamodel.readJSONheadedASCII(str(hr_path))
            # print(f'file_name={hr_path}, L-shell_keys={[key for key in hr.keys() if "wainL" in key]}')
            hr['Time'] = pd.to_datetime(hr['Time'])
            try:
                cadence = float(hr.attrs['CADENCE'])
            except KeyError as err:
                if 'CADENCE' in str(err):
                    print(f'hr_path={hr_path}')
                    cadence = input(f'What is the FU{self.sc_id} cadence for {hr_path.name}?')
                    cadence = float(cadence)
                else:
                    raise
                
            # All of the code to detect microbursts is here.
            s = signal_to_background.FirebirdSignalToBackground(
                hr['Col_counts'], cadence, 
                self.background_width_s, 
                self.microburst_width_s
                )
            s.significance()
            try:
                s.find_microburst_peaks(std_thresh=self.std_thresh, 
                                        detect_channel=self.detect_channel)
            except ValueError as err:
                if str(err) == 'No detections found':
                    continue
                else:
                    raise
                
            daily_microburst_list = pd.DataFrame(
                data=np.nan*np.ones((len(s.peak_idt), len(self.catalog_columns)), dtype=object), 
                columns=self.catalog_columns
                )
            daily_microburst_list.loc[:, self.hr_keys] = np.array(
                [hr[col][s.peak_idt] for col in self.hr_keys],
                dtype=object
                ).T
            daily_microburst_list.loc[:, self.count_keys] = hr['Col_counts'][s.peak_idt, :]
            daily_microburst_list.loc[:, self.sig_keys] = s.n_std.loc[s.peak_idt, :].to_numpy()
                                            
            self.microburst_list = pd.concat((self.microburst_list, daily_microburst_list))
            
        self.microburst_list = self.microburst_list.reset_index()
        del(self.microburst_list['index'])
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
            

if __name__ == '__main__':
    sc_id = 4
    microburst_width_s = 0.1
    background_width_s = 0.5
    std_thresh = 10

    s = SignalToBackgroundLoop(sc_id, microburst_width_s, background_width_s, std_thresh)
    s.loop()
    s.save_microbursts()