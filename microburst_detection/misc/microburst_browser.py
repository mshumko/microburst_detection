import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import pandas as pd
from datetime import date, datetime, timedelta
import dateutil.parser
from matplotlib.widgets import Button, TextBox
import pathlib

from microburst_detection.misc.load_firebird import readJSONheadedASCII
from microburst_detection import config

class Browser:
    def __init__(self, fb_id, 
                catalog_name=None, plot_width=5, 
                catalog_save_name=None, filterDict={}, 
                jump_to_latest=True):
        """
        This class plots the FIREBIRD microbursts and allows the user to browse
        detections in the future and past with buttons. Also there is a button
        to mark the event as a microburst.

        DOES NOT WORK WITH NEWER VERSIONS OF MATPLOTLIB DUE TO A RACE CONDITION
        Issue: https://github.com/matplotlib/matplotlib/issues/9491/
        Fix: works with matplotlib-2.2.3 --  set up pyenv to maintain different
        versions.
        """
        self.fb_id = fb_id
        self.plot_width = plot_width

        self.load_catalog(catalog_name=catalog_name)
        self.filter_catalog(filterDict=filterDict)

        if catalog_save_name is None:
            self.catalog_save_name = (f'{catalog_name.split(".")[0]}_sorted.txt')
        else:
            self.catalog_save_name = catalog_save_name
        self.catalog_save_path = pathlib.Path(
            config.PROJECT_DIR, 
            'data',                     
            self.catalog_save_name)

        # Load the filtered catalog if it already exists. This is
        # userful if you can't sort all microbursts at once!
        if os.path.exists(self.catalog_save_path):
            self.load_filtered_catalog_indicies()
        else:
            self.microburst_idx = np.array([])

        self.current_date = date.min
        self._init_plot()
        if jump_to_latest and len(self.microburst_idx):
            self.index = self.microburst_idx[-1]
        else:
            # Start at row 0 in the dataframe.
            self.index = 0 
        self.plot()
        return

    def load_catalog(self, catalog_name):
        """
        Load the catalog and convert the times
        """
        catalog_path = pathlib.Path(
            config.PROJECT_DIR, 
            'data', catalog_name
            )
        self.catalog = pd.read_csv(catalog_path)
        self.catalog.Time = pd.to_datetime(self.catalog.Time)
        return

    def filter_catalog(self, filterDict={}):
        """ Filter the catalog with a default filter + keys in filterDict """
        for key in filterDict:
            if hasattr(filterDict[key], '__len__'):
                self.catalog = self.catalog[(self.catalog[key] >= filterDict[key][0]) &
                                             (self.catalog[key] <= filterDict[key][1])]
            else:
                self.catalog = self.catalog[(self.catalog[key] == filterDict[key])]
        return

    def next(self, event):
        """ Plots the next detection """
        # Just return if at the end of the dataframe.
        if self.index + 1 >= self.catalog.shape[0]:
            return
        self.index += 1
        self.plot()
        return

    def prev(self, event):
        """ Plots the previous detection """
        # Just return if at the end of the dataframe.
        if self.index == 0:
            return
        self.index -= 1
        self.plot()
        return

    def append_remove_microburst(self, event):
        """ 
        Appends or removes the current catalog row to 
        self.filtered_catalog which will then
        be saved to a file for later processing.
        """
        if self.index not in self.microburst_idx:
            self.microburst_idx = np.append(self.microburst_idx, self.index)
            self.bmicroburst.color = 'g'
            print('Microburst saved at', self.catalog.iloc[self.index].Time)
        else:
            self.microburst_idx = np.delete(self.microburst_idx, 
                np.where(self.microburst_idx == self.index)[0])
            self.bmicroburst.color = '0.85'
            print('Microburst removed at', self.catalog.iloc[self.index].Time)
        return

    def key_press(self, event):
        """
        Calls an appropriate method depending on what key was pressed.
        """
        if event.key == 'm':
            # Mark as a microburst
            self.append_remove_microburst(event)
        elif event.key == 'a':
            # Move self.index back and replot.
            self.prev(event)
        elif event.key =='d':
            # Move the self.index forward and replot.
            self.next(event)
        return
       
    def change_index(self, index):
        self.index = int(index)
        self.plot()
        return

    def plot(self):
        """ 
        Given a self.current_row in the dataframe, make a space-time plot 
        """
        print('Index position = {}/{}'.format(
                    self.index, self.catalog.shape[0]-1))
        current_row = self.catalog.iloc[self.index]
        self.index_box.set_val(self.index)
        #self._clear_ax()
        self.ax.clear()

        if current_row['Time'].date() != self.current_date:
            # Load current day data if not loaded already
            print('Loading data from {}...'.format(current_row['Time'].date()), 
                    end=' ', flush=True)
            self._load_hr(current_row['Time'].date())
            self.current_date = current_row['Time'].date()
            print('done.')

        # Turn microburst button green if this index has been marked as a microburst.
        if self.index in self.microburst_idx:
            self.bmicroburst.color = 'g'
        else:
            self.bmicroburst.color = '0.85'

        dt = timedelta(seconds=self.plot_width/2)
        time_range = (current_row['Time'] - dt, current_row['Time'] + dt)
        idt = np.where((self.hr['Time'] > time_range[0]) & 
                        (self.hr['Time'] < time_range[1]))[0]
        #self.make_plot(current_row, savefig=False)
        self.ax.plot(self.hr['Time'][idt], self.hr['Col_counts'][idt, :])
        self.ax.axvline(current_row['Time'])
        self.ax.set_title('FIREBIRD Microburst Browser\n {} {}'.format(
                        current_row['Time'].date(), 
                        current_row['Time'].time()))
        self.ax.set_ylabel(f'Col counts/{self.cadence} ms')
        self.ax.set_xlabel('UTC')
        # if isinstance(current_row, pd.DataFrame):
        self._print_aux_info(current_row)
        self.ax.set_xlim(*time_range)
        plt.draw()
        return

    def _print_aux_info(self, current_row):
        """ Print separation info as well as peak width info to the canvas. """
        self.textbox.clear()
        self.textbox.axis('off')
        col1 = (f'L = {round(current_row.McIlwainL, 1)}\n'
                f'MLT = {round(current_row.MLT, 1)}\n'
                f'Lat = {round(current_row.Lat, 1)}\n'
                f'Lon = {round(current_row.Lon, 1)}\n')  
        self.textbox.text(0, 1, col1, va='top')
        #self.textbox.text(1.3, 1, col2, va='top')
        return

    def _clear_ax(self):
        [a.clear() for a in self.ax]
        return 

    def _init_plot(self):
        """
        Initialize subplot objects and text box.
        """
        fig, self.ax = plt.subplots(1, figsize=(8, 7))
        plt.subplots_adjust(bottom=0.2)

        # Define button axes.
        self.axprev = plt.axes([0.54, 0.06, 0.12, 0.075])
        self.axburst = plt.axes([0.67, 0.06, 0.13, 0.075])
        self.axnext = plt.axes([0.81, 0.06, 0.12, 0.075])

        # Define buttons and their actions.
        self.bnext = Button(self.axnext, 'Next (d)', hovercolor='g')
        self.bnext.on_clicked(self.next)
        self.bprev = Button(self.axprev, 'Previous (a)', hovercolor='g')
        self.bprev.on_clicked(self.prev)
        self.bmicroburst = Button(self.axburst, 'Microburst (m)', hovercolor='g')
        self.bmicroburst.on_clicked(self.append_remove_microburst)

        # Define the textbox axes.
        self.textbox = plt.axes([0.1, 0.05, 0.2, 0.075])
        self.textbox.axis('off')
        # Define index box.
        self.axIdx = plt.axes([0.59, 0.01, 0.32, 0.04])
        self.index_box = TextBox(self.axIdx, 'Index')
        self.index_box.on_submit(self.change_index)

        # Initialise button press
        fig.canvas.mpl_connect('key_press_event', self.key_press)
        return

    def save_filtered_catalog(self):
        """
        For every index that a user clicked microburst on, save
        those rows from the catalog into a new catalog with the
        name of self.catalog_save_name.
        """
        # Return if there are no micriobursts to save.
        if not hasattr(self, 'microburst_idx'):
            return
        # Remove duplicates indicies
        self.microburst_idx = np.unique(self.microburst_idx)
        #save_path = os.path.join(catalog_save_dir, self.catalog_save_name)
        print('Saving filtered catalog to {}'.format(self.catalog_save_path))
        df = self.catalog.iloc[self.microburst_idx]
        df.to_csv(self.catalog_save_path, index=False)
        return

    def load_filtered_catalog_indicies(self):
        """
        Load a filtered catalog and populate the self.microbirst_idx array
        with existing detections. This method exists to help the user resume
        the 
        """
        filtered_catalog = pd.read_csv(self.catalog_save_path)
        filtered_catalog.Time = pd.to_datetime(filtered_catalog.Time)
        catalog_times_numeric = date2num(self.catalog.Time.dt.to_pydatetime())
        filtered_catalog_times_numeric = date2num(filtered_catalog.Time.dt.to_pydatetime())
        self.microburst_idx = np.where(np.in1d(catalog_times_numeric, 
                                    filtered_catalog_times_numeric, 
                                    assume_unique=True))[0]
        return

    def _load_hr(self, date):
        """ Loads FIREBIRD HiRes data and converts times """
        search_str = f'FU{self.fb_id}_Hires_{date}_L2.txt'
        hr_paths = list(pathlib.Path(config.FB_DIR).rglob(search_str))
        assert len(hr_paths) == 1, (f'A unique HiRes path not found.\n'
                                    f'hr_paths={hr_paths} search_str={search_str}')
        hr_path = str(hr_paths[0])
        self.hr = readJSONheadedASCII(hr_path)
        self.hr['Time'] = pd.to_datetime(self.hr['Time'])
        #self.hr['Time'] = np.array([dateutil.parser.parse(t) for t in self.hr['Time']])
        self.cadence = 1000*float(self.hr.attrs['CADENCE'])
        return

### ARLO'S CODE ###
sc_id = 4
callback = Browser(sc_id, catalog_name=f'FU{sc_id}_microbursts.csv')

# Initialize the GUI
plt.show()
# Save the catalog.
callback.save_filtered_catalog()
