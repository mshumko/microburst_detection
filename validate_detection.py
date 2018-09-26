# This code will explore the microburst detecting properties of the burst 
# parameter and wavelet methods.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import burst_parameter.burst_parameter as burst_parameter
#import waveletAnalysis

def gaus(t, *args):
    """
    Returns the value of a Gaussian of the form 
    f(t, args) = args[0]*exp(-(t - args[1])^2/(2*args[2]^2))
    Multiple of 3 parameters in args can be given. For each triplet,
    the first arg is the amplitude, second is the t_0, and third is
    sigma.
    """
    args = args[0]
    assert len(args) % 3 == 0, 'Length of arguments must be a factor of 3!'
    val = 0
    for i in range(0,len(args),3): # Loop over arg triplets.
        x = np.divide(np.power((t-args[i+1]), 2),(2*np.power(args[i+2], 2)))
        val += args[i]*np.exp(-x.astype(float))
    return val


if __name__ == '__main__':
    #nTrials = 5
    dT = 0.1
    t = np.arange(-5, 5, dT)
    cadence = 0.1
    t0 = 0
    A = 1000
    baseline = np.linspace(0, 1000)
    fwhm = 0.1