import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import mc_model_config

def gaus(x, p):
    """ 
    A plain old Gaussian function. The amplitude, 
    offset, and sigma are wrapped up into the p 
    tuple argument.
    """
    A, x0, sigma = p
    return A*np.exp(-0.5*((x-x0)/sigma)**2)

class ModelPeaks:
    def __init__(self, config):
        """
        A Monte Carlo (MC) model that generates a bunch of Gaussian profiles
        with Poisson noise.
        """
        self.config = config
        self.n_data_points = int(self.config.time_width_s/self.config.cadence_s)
        self.time_array = np.arange(0, self.config.time_width_s, self.config.cadence_s)

        # Run the MC count simulation
        self.simulate_peaks()
        return

    def simulate_peaks(self, add_noise=True):
        """
        Generate n_iter number of Gaussian profiles with 
        Poisson noise.
        """
        # Make array of n_data_points x n_iter
        # filled with the background amplitude counts.
        ones_array = np.ones((self.n_data_points, self.config.n_iter), dtype=int)
        self.counts = self.config.background_a*ones_array
        # Generate an array of peak widths (half the width ~ std, will be fed into gaus)
        self.peak_widths = self.config.peak_width_dist.rvs(size=self.config.n_iter)

        # Generate n_iter number of Gaussians and add to background 
        for i in range(self.config.n_iter):
            p = (
                self.config.peak_a, 
                self.config.time_width_s/2, 
                self.peak_widths[i]/2 # Half the peak width to get ~std
                )
            self.counts[:, i] += np.array(gaus(self.time_array, p)).astype(int)

        # Add Poisson noise
        if add_noise:
            self.counts = np.random.poisson(lam=self.counts)
        return

    def visualize_model(self, n_plot=10, ax=None):
        """
        Visualize a subset of the random Gaussian profiles.
        """
        idx_plot = np.random.randint(0, self.counts.shape[0], size=n_plot)

        if ax is None:
            _, ax = plt.subplots(2)
        
        for counts_i in self.counts[:, idx_plot].T:
            ax[0].plot(self.time_array, counts_i, 'k')

        # Histogram the peak widths
        ax[1].hist(self.peak_widths, histtype='step', lw=3, color='k')
        
        ax[0].set(xlabel='Time', ylabel='Counts', 
            title='Microburst Width Monte Carlo Model')
        ax[1].set(xlabel='Peak width', ylabel='Number of peaks',
            title='Distribtuion of peak widths')
        plt.tight_layout()
        return

if __name__ == "__main__":
    m = ModelPeaks(mc_model_config)
    m.visualize_model()
    plt.show()