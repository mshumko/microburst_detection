# My own version of spacepy.datamodel.readJSONheadedASCII
import json

import numpy as np
import pandas as pd


def readJSONheadedASCII(file_path):
    """
    My simple implementation of spacepy.datamodel.readJSONheadedASCII that
    is specific for FIREBIRD-II data. You may use this if you can't install 
    spacepy for whatever reason.
    """
    # Read in the JSON header.
    header_list = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_list.append(line[1:])
            else:
                raw_data_str = f.readlines()

    # Massage the header
    clean_header_str = ''.join(header_list).replace('\n', '')
    parsed_header = json.loads(clean_header_str)
    # Massage the data
    raw_data_str = [row.replace('\n', '') for row in raw_data_str]
    # Parse times
    times_str = [row.split()[0] for row in raw_data_str]
    # Parse the other columns
    data_converted = np.array([row.split()[1:] for row in raw_data_str]).astype(float)

    data = HiRes()
    data['Time'] = pd.to_datetime(times_str)
    for key in parsed_header:
        key_header = parsed_header[key]
        data.attrs[key] = key_header  # Save the attribute data.

        if key == 'Time':  # We already added Time.
            continue
        # Header key that correspond to columns
        if isinstance(key_header, dict):
            if len(key_header['DIMENSION']) != 1:
                raise NotImplementedError(
                    "readJSONheadedASCII doesn't implement columns with more than "
                    f"1 multidimensional. Got {key_header['DIMENSION']}."
                    )
            start_column = key_header['START_COLUMN']-1
            end_column = key_header['START_COLUMN']-1+key_header['DIMENSION'][0]
            if key_header['DIMENSION'][0] == 1:
                data[key] = data_converted[:, start_column]
            else:
                data[key] = data_converted[:, start_column:end_column]
        else:
            # Header key that correspond to global attributes
            if key in ['CADENCE', 'CAMPAIGN']:
                data.attrs[key] = float(key_header)
            else:
                data.attrs[key] = key_header
    return data

class HiRes(dict):
    """
    Expand Python's dict class to include an attr attribute dictionary.
    
    Code credit goes to Matt Anderson:
    https://stackoverflow.com/questions/2390827/how-to-properly-subclass-dict-and-override-getitem-setitem
    (blame him for problems)
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        self.attrs = {}
        return
