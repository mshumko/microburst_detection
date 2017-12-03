# This code will explore the microburst detecting properties of the burst 
# parameter and wavelet methods.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
        colors = cm.rainbow(np.linspace(0, 1, num=self.ts.shape[1]))

        if ax is None:
            fig, self.ax = plt.subplots()
        else:
            self.ax = ax

        for (i, color) in enumerate(colors):
            self.ax.plot(self.t, self.ts[:, i], c=color)

        self.ax.set(title='Random Gaussians', xlabel='Time (s)', ylabel='Counts/s')

        if ax is None:
            plt.show()
        return


class TestBurstParam(CreateFakeData):
    def __init__(self, A, t, t0):

        return





if __name__ == '__main__':
    nTrials = 10
    dT = 0.1
    t = np.arange(-5, 5, dT)
    fwhm = np.linspace(0.1, 2, num=nTrials)
    t0 = 0
    A = 100
    fwhm = np.linspace(0.1, 2, num=nTrials)
    fakeData = CreateFakeData(t)
    fakeData.makeTimeseries(A, t0, fwhm, noise=True)
    fakeData.plotTimeseires()