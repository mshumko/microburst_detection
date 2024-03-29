# Introduction
This repo contains a few algorithms to detect electron microbursts. The two 
implemented methods are:
- Signal to background
- Wavelet-based filtering 

The two methods involve some sort of count or standard deviation threshold
and must be tuned for your data. The signal to background method is 
mathematically the same as Paul O'Brien's burst parameter but with 
different parameters The burst parameter was used in, for example, 
O'Brien et al., 2003; Douma et al., 2017; Shumko et al., 2020.

# Dependencies:
Easiest method to install the numpy, pandas, and progressbar 
dependencies is with pip and the requirements.txt file. 

Install with: 
```python3 -m pip install -r requirements.txt```

And set the project and data paths generally with:
```python3 -m microburst_detection init```

Otherwise if you're a FIREBIRD user:
```python3 -m microburst_detection firebird```

# Example Uses:

## Signal to background detection
To play around with this detection method the script in the ```python3 if __name__ == '__main__'``` block of ```signal_to_background.py``` has the parameters to tweak this detector, change the data analyzed, and quickly visualize the results. The two classes are very similar, one works with generic 1d count time series (from one energy channel), while the other works with 2d count time series where the two dimensions are: nTime x nEnergyChannels.

![Example of the signal-to-background program flagging microbursts](/example_plots/signal_to_background_example.png)

The ```signal_to_background_loop.py``` calls ```signal_to_background.py``` on all FIREBIRD HiRes data and saves it to a csv file in ```<<project_folder>>/data/``` folder where ```<<project_folder>>``` is specified in dirs.py.

## Wavelet-bassed Microburst Detection
The other microburst detection method is based on wavelet filtering in the frequency-time domain. This method is heavily based on the [Torrence and Compo, 1998](https://psl.noaa.gov/people/gilbert.p.compo/Torrence_compo1998.pdf) paper and the wavelet analysis code is adapted from their [GitHub repo](https://github.com/chris-torrence/wavelets)

![Wavelet microburst detection](/example_plots/wavelet_detection_example.png)

In a nutshell the wavelet-based method includes the following steps.
1. Transform the time series into the wavelet domain using a predetermined wavelet basis, as shown in Panels (a) -> (b). (The Mexican Hat wavelet basis by default).
2. Filter out the statistically insignificant wavelet power using a red-noise spectrum. 
3. Filter the wavelet power by oscillation period (panel (b) shows the white-hached power spectrum that is filtered out).
4. Inverse wavelet transform the remaining wavelet spectrum and you get the filtered counts, in arbitrary units in Panel (c).
5. Apply a threshold test. Anything above the blue dashed line threshold in Panel C is a candidate microburst. Typically more than one data point is above the threshold for each interval, so every continuous interval above the threshold is identified by the ```locate_consecutive_numbers``` program and in each interval the peak count rate is found. 

More details are provided by [Torrence and Compo, 1998](https://psl.noaa.gov/people/gilbert.p.compo/Torrence_compo1998.pdf)


## Manually sort microbursts
Once you ran ```signal_to_background_loop.py``` and the catalog file is generated, you can use the microburst browser located in ```misc/microburst_browser.py``` GUI to sort the microburst detections. In this GUI you can navigate forward and backward in the list and mark microbursts. Besides the navigation buttons you can also use your keyboard keys to navigate. If you accidently mark something as microburst, press it again and it will be removed. 

![Microburst browser](/example_plots/microburst_browser.png)

Upon exiting the browser the marked microbursts will be saved to where the original catalog was located, suffixed with "_sorted". If you open the browser again and the sorted catalog is populated, the browser will jump to the most recent microburst detection. Therefore you can stop and pick up your sorting from the last microburst detection.

The last feature lets you jump to any microburst in the list by changing the index in the GUI's index box in the bottom right.

# Bibliography

Blum, L., X. Li, and M. Denton (2015), Rapid MeV electron precipitation as observed by SAMPEX/HILT during high-speed stream-driven storms, J. Geophys. Res. Space Physics, 120, 3783–3794, doi:10.1002/2014JA020633.

Douma, E., C. J. Rodger, L. W. Blum, and M. A. Clilverd (2017), Occurrence characteristics of relativistic electron microbursts from SAMPEX observations, J. Geophys. Res. Space Physics, 122, doi:10.1002/2017JA024067.

Greeley, A. D., Kanekal, S. G., Baker, D. N., Klecker, B., & Schiller, Q. (2019). Quantifying the contribution of microbursts to global electron loss in the radiation belts. Journal of Geophysical Research: Space Physics, 124. https://doi.org/10.1029/2018JA026368

O’Brien, T. P., K. R. Lorentzen, I. R. Mann, N. P. Meredith, J. B. Blake, J. F. Fennell, M. D. Looper, D. K. Milling, and R. R. Anderson, Energization of relativistic electrons in the presence of ULF power and MeV microbursts: Evidence for dual ULF and VLF acceleration, J. Geophys. Res., 108(A8), 1329, doi:10.1029/2002JA009784, 2003.

Shumko, M., Johnson, A. T., Sample, J. G., Griffith, B. A., Turner, D. L., O'Brien, T. P., et al. (2020). Electron microburst size distribution derived with AeroCube-6. Journal of Geophysical Research: Space Physics, 125, e2019JA027651. https://doi.org/10.1029/2019JA027651