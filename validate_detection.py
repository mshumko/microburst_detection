# This code will explore the microburst detecting properties of the burst 
# parameter and wavelet methods.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import burst_parameter.burst_parameter as burst_parameter
#import waveletAnalysis

class CreateFakeData:
    def __init__(self, t):
        self.t = t # Time array (not datetimes)
        return

    @staticmethod
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

    def makeTimeseries(self, A, t0, widthArr, noise=True, baseline=1):
        """
        This function creates an nTrials numer of timeseires in an
        nT x len(widthArr) sized array that contains a Gaussian of amplitude
        A and center t0. The widthArr is an array that contains the FWHM of
        each Gaussian.

        If noise is true, Poisson noise will be added.
        The timeseries is shifted up by baseline variable.
        """
        self.A = A
        self.baseline = baseline

        sigmaArr = widthArr/2.3548 # Convert FWHM to sigmas.
        self.ts = baseline*np.ones((len(self.t), len(widthArr)), dtype=float)

        for (i, s) in enumerate(sigmaArr):
            self.ts[:, i] += self.gaus(self.t, (A, t0, s))

        if noise:
            self.ts = np.random.poisson(lam=self.ts)
        return

    def plotTimeseires(self, nth=1, ax=None):
        """
        This function will plot the Gaussian timeseries. nth=1 is an optional
        parameter that is used to downsample the self.ts array plot every nth 
        Gaussian.
        """
        if ax is None:
            fig, self.ax = plt.subplots()
        else:
            self.ax = ax

        colors = cm.rainbow(
            np.linspace(0, 1, num=self.ts.shape[1])) # For plotting.
        for (i, color) in enumerate(colors):
            self.ax.plot(self.t, self.ts[:, i], c=color)

        self.ax.set(title='Random Gaussians | A = {} | baseline = {}'.format(
            self.A, self.baseline), xlabel='Time (s)', ylabel='Counts/s')

        if ax is None:
            plt.show()
        return


class TestBurstParam(CreateFakeData):
    def __init__(self, A, t, t0, widthArr, baseline=1):
        CreateFakeData.__init__(self, t)
        self.makeTimeseries(A, t0, widthArr, baseline=baseline)
        return

    def calcParam(self, cadence, N_WIDTH=0.1, A_WIDTH=0.5):
        self.detectionParam = np.nan*np.ones_like(self.ts)

        for (i, counts) in enumerate(self.ts):
            self.detectionParam[i, :] = burst_parameter.obrien_burst_param(
                counts, cadence, N_WIDTH=N_WIDTH, A_WIDTH=A_WIDTH)
        return

    def calcThresh(self, thresh):
        """ 
        This function will calculate the burst parameter points that are above the 
        thresh parameter.
        """
        self.thresh = thresh
        self.detectInd = np.where(self.detectionParam > thresh)
        return

    def plotDataDetection(self, ax=None):
        """
        This function will plot the Gaussian timeseries. nth=1 is an optional
        parameter that is used to downsample the self.ts array plot every nth 
        Gaussian.
        """
        if ax is None:
            fig, self.bx = plt.subplots(figsize=(11, 8))
        else:
            self.bx = ax

        colors = cm.rainbow(
            np.linspace(0, 1, num=self.ts.shape[1])) # For plotting.
        for (i, color) in enumerate(colors): # Plot burst parameter
            self.bx.plot(self.t, self.detectionParam[:, i], c=color)

        # Plot detections
        if hasattr(self, "detectInd"):
            tt = np.repeat(self.t[:, np.newaxis], self.ts.shape[1], axis=1)
            self.bx.scatter(tt[self.detectInd], 
                self.detectionParam[self.detectInd], c='k')
            self.bx.axhline(self.thresh)
            self.bx.text(self.t[0], self.thresh, 'Threshold', va='bottom',
                )

        self.bx.set(ylabel='Burst Parameter', xlabel='Time (s)')

        if ax is None:
            plt.show()
        return


if __name__ == '__main__':
    nTrials = 10
    dT = 0.1
    t = np.arange(-5, 5, dT)
    cadence = 0.1
    t0 = 0
    A = 10
    baseline=1
    fwhm = np.linspace(0.1, 2, num=nTrials)

    testObj = TestBurstParam(A, t, t0, fwhm, baseline=baseline)
    testObj.calcParam(cadence)
    testObj.calcThresh(2)
    
    ### Plotting ###
    fig, ax = plt.subplots(2, sharex=True)
    testObj.plotTimeseires(ax=ax[0])
    testObj.plotDataDetection(ax=ax[1])
    plt.tight_layout()
    plt.show()