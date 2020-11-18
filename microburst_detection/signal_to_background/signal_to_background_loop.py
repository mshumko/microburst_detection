import numpy as np
import pandas as pd
import pathlib
import progressbar

import spacepy

import signal_to_background
import microburst_detection.dirs

class SignalToBackgroundLoop:
    def __init__(self, sc_id, background_width_s, std_thresh, detect_channel=0):
        """
        This program uses signal_to_background detection code to
        loop over all of the FIREBIRD data and detect all 
        microbursts.

        Parameters
        ----------
        sc_id : int
            Spacecraft id. Either 3 or 4
        background_width_s : float
            The baseline width in time to calculate the running mean
        std_thresh : float
            The baseline standard deviation threshold above the baseline
            that the data point must be to satisfy the microburst criteria
        detect_channel : int, optional
            The FIREBIRD energy channel number to use for std_thresh 
            critera. This is channel 0 by default.
        """
        self.sc_id = sc_id
        self.background_width_s = background_width_s
        self.std_thresh = std_thresh
        self.detect_channel = detect_channel

        # Find all of the HiRes files
        search_str = f'FU{sc_id}_Hires_*.txt'
        hr_dir = microburst_detection.dirs.firebird_dir(sc_id)
        self.hr_paths = sorted(pathlib.Path(hr_dir).rglob(search_str))
        return

    def loop(self, save_keys='default'):
        """
        Loop over all the HiRes data and run the signal_to_background
        microburst detector on every day. For the detected microbursts
        save a handful of columns specified by the save_keys kwarg to
        self.microburst_list.

        If save_keys='default', a default set of keys will be used
        including time, the spatial, and geomagnetic coordinates. These 
        keys must be in the HiRes data. Furthermore if you're running
        this on the FIREBIRD data, the baseline standard deviation values
        for the 6 energy channels will be saved as well.
        """
        if save_keys == 'default':
            save_keys = (['Time', 'Lat', 'Lon', 'Alt', 
                        'McIlwainL', 'MLT', 'kp'] + 
                        [f'sig{i}' for i in range(6)])
        # hr_keys only has the keys that are in the HiRes data.
        hr_keys = [col for col in save_keys if 'sig' not in col]
        
        #self.microburst_list = pd.DataFrame(columns=save_keys)
        self.microburst_list = np.nan*np.ones(
                (0, len(save_keys)), dtype=object
                )

        for hr_path in progressbar.progressbar(self.hr_paths):
            hr = spacepy.datamodel.readJSONheadedASCII(str(hr_path))
            hr['Time'] = pd.to_datetime(hr['Time'])
            cadence = float(hr.attrs['CADENCE'])

            # All of the code to detect microbursts is here.
            s = signal_to_background.FirebirdSignalToBackground(
                hr['Col_counts'], cadence, 
                self.background_width_s
                )
            s.significance()
            try:
                s.find_microburst_peaks(std_thresh=self.std_thresh, 
                                        detect_channel=detect_channel)
            except ValueError as err:
                if str(err) == 'No detections found':
                    continue
                else:
                    raise

            daily_detections = np.nan*np.ones(
                (len(s.peak_idt), len(save_keys)), dtype=object
                )
            # Loop over each microburst detection and append to 
            # the daily list
            for i, peak_i in enumerate(s.peak_idt):
                # Save certain HiRes heys
                #print(s.n_std.shape)
                daily_detections[i, :len(hr_keys)] = [hr[col][peak_i] for col in hr_keys] 
                # Save the significance above baseline values for all channels.
                daily_detections[i, len(hr_keys):] = s.n_std.loc[peak_i, :]
                                            
            self.microburst_list = np.concatenate((self.microburst_list, daily_detections))

        self.microburst_list = pd.DataFrame(data=self.microburst_list, 
                                            columns=save_keys)
        return

    def save_microbursts(self, save_name=None):
        """
        Save the microburst list to a csv file in the data/ directory.
        If the directory does not exist, one will be created.
        """
        save_dir = pathlib.Path(
            microburst_detection.dirs.project_dir,
            'data'
        )
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
    background_width_s = 2
    std_thresh = 10

    s = SignalToBackgroundLoop(sc_id, background_width_s, std_thresh)
    s.loop()
    s.save_microbursts()
