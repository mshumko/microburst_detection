import matplotlib.pyplot as plt
import numpy as np
import progressbar

import signal_to_background.signal_to_background as signal_to_background
import mc_model
import mc_model_config

m = mc_model.ModelPeaks(mc_model_config) # Generate a bunch of time series.

detected = np.ones(m.peak_widths.shape[0], dtype=bool)

background_width_s = 10
sig_thresh_std = 10

for i in progressbar.progressbar(range(m.counts.shape[1])):
    counts = m.counts[:, i]
    s = signal_to_background.SignalToBackground(
        counts, mc_model_config.cadence_s, 
        background_width_s
    )
    s.significance()

    detected_points = np.where(s.n_std >= sig_thresh_std)[0]
    if len(detected_points):
        detected[i] = True
    else:
        detected[i] = False


def visualize_peaks(time, counts, ax, n_plot=20):
    """
    Visualize a subset of the random Gaussian profiles.
    """
    idx_plot = np.random.randint(0, counts.shape[1], size=n_plot)
    
    for counts_i in counts[:, idx_plot].T:
        ax.plot(time, counts_i, 'k')
    return

fig, ax = plt.subplots(3, figsize=(6, 7), sharex=True)

fig2, bx = plt.subplots()

bins = np.arange(0, 11)
ax[0].hist(m.peak_widths, bins=bins, histtype='step', lw=3, color='k')
ax[1].hist(m.peak_widths[detected], bins=bins, histtype='step', lw=3, color='k')
ax[2].hist(m.peak_widths[~detected], bins=bins, histtype='step', lw=3, color='k')

y_lims = ax[0].get_ylim()

ax[0].set(ylabel='True number of peaks', 
        title='Signal-to-baseline Detector Sensitivity\nMonte Carlo Model',
        ylim=(None, 1.5*y_lims[1]))
ax[1].set(ylabel='# of peaks detected')
ax[2].set(xlabel='Peak width [s]', ylabel='# of peaks not detected')

annotate_str = (f'background_width = {background_width_s} s\n'
                f'std_thresh = {sig_thresh_std}\n'
                f'peak_amplitude = {mc_model_config.peak_a}\n'
                f'background_amplitude = {mc_model_config.background_a}\n')
ax[0].text(0.99, 0.99, annotate_str, transform=ax[0].transAxes, ha='right', va='top')

visualize_peaks(m.time_array, m.counts, bx)
bx.set(xlabel='Time', ylabel='Counts')


plt.tight_layout()
plt.show()