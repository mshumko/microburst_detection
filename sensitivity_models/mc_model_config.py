import scipy.stats

### Instrument settings
# Total time width of the time sereis (needs to fully encompass the detection analysis) 
time_width_s = 60
# Instrument cadence
cadence_s = 0.1

### Peak settings
# background amplitude
background_a = 0 
# Peak amplitude
peak_a = 100 
# uniform peak widths between 0 and 10 seconds.
peak_width_dist = scipy.stats.uniform(loc=0, scale=5) 

### Simulation settings
#n_widths = 100 # Take n_widths from peak_width_dist
# Generate n_iter time series with noise, with a peak width 
# randomly picked from peak_width_dist
n_iter = 10000
