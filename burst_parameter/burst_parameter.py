import numpy as np

def moving_average(interval, window_size):
    """
    NAME: moving_average(interval, window_size)
    USE: Calculates a running average of size window_size on the interval array.
    INPUTS: interval: an 1d array of data.
            window_size: amount of points to average over.
    OUTPUTS: Average of the interval array over a window_size length.
    AUTHOR: Mykhaylo Shumko
    MOD: 2017-10-28
    """
    window= np.ones(window_size)/float(window_size)
    avg = np.convolve(interval, window, 'same')
    return avg

def obrien_burst_param(counts, cadence, N_WIDTH = 0.100, A_WIDTH = 0.500):
    """
    NAME: obrien_burst_param(counts, cadence, N_WIDTH = 0.100, A_WIDTH = 0.500)
    USE: Calculates Paul O'brien's microburst parameter on count data. Refer to
         O'Brien et al. 2003 JGR.
    INPUTS: REQUIRED:
            counts: an 1d array of count data.
            cadence: Data cadence in the same units as n_width and a_width.
            OPTIONAL:
            n_width = 0.100 : Time width of the n parameter, refered to as n100.
            a_width = 0.500 : Time width of the a parameter, refered to as a500.
    OUTPUTS: Burst parameter array of the same size as counts..
    AUTHOR: Mykhaylo Shumko
    MOD: 2017-10-28
    """
    n = movingaverage(counts, N_WIDTH//cadence)
    a = movingaverage(counts, A_WIDTH//cadence)
    bp = (n-a)/np.sqrt(1+a)
    return bp
