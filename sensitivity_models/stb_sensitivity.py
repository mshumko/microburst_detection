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

fig, ax = plt.subplots(3, sharex=True, figsize=(6, 10))

bins = np.arange(0, 11)
ax[0].hist(m.peak_widths, bins=bins, histtype='step', lw=3, color='k')
ax[1].hist(m.peak_widths[detected], bins=bins, histtype='step', lw=3, color='k')
ax[2].hist(m.peak_widths[~detected], bins=bins, histtype='step', lw=3, color='k')

ax[0].set(ylabel='True number of peaks', 
        title='Signal-to-baseline Sensitivity\nMonte Carlo Model')
ax[1].set(ylabel='# of peaks detected')
ax[2].set(xlabel='Peak width [s]', ylabel='# of peaks not detected')

plt.show()