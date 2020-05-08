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
Easiest method to install the numpy, pandas, progressbar, and spacepy 
dependencies is with pip and the requirements.txt file. 

Example: ```pip3 install -r requirements.txt```

# Example Uses:


# Project Structure
├── burst_parameter
│   └── burst_parameter.py
├── data
│   └── FU4_microbursts.csv
├── detect_microbursts.py
├── dirs.py
├── example_plots
├── misc
│   ├── locate_consecutive_numbers.py
│   └── microburst_browser.py
├── README.md
├── requirements.txt
├── sensitivity_models
│   ├── mc_model_config.py
│   ├── mc_model.py
│   └── stb_sensitivity.py
├── signal_to_background
│   ├── signal_to_background_loop.py
│   └── signal_to_background.py
├── validate_detection.py
├── validate_sampex.py
└── wavelets
    ├── waveletAnalysis.py
    ├── waveletFunctions.py
    └── wavelet_test_scipt.p

# Bibliography

Blum, L., X. Li, and M. Denton (2015), Rapid MeV electron precipitation as observed by SAMPEX/HILT during high-speed stream-driven storms, J. Geophys. Res. Space Physics, 120, 3783–3794, doi:10.1002/2014JA020633.

Douma, E., C. J. Rodger, L. W. Blum, and M. A. Clilverd (2017), Occurrence characteristics of relativistic electron microbursts from SAMPEX observations, J. Geophys. Res. Space Physics, 122, doi:10.1002/2017JA024067.

Greeley, A. D., Kanekal, S. G., Baker, D. N., Klecker, B., & Schiller, Q. (2019). Quantifying the contribution of microbursts to global electron loss in the radiation belts. Journal of Geophysical Research: Space Physics, 124. https://doi.org/10.1029/2018JA026368

O’Brien, T. P., K. R. Lorentzen, I. R. Mann, N. P. Meredith, J. B. Blake, J. F. Fennell, M. D. Looper, D. K. Milling, and R. R. Anderson, Energization of relativistic electrons in the presence of ULF power and MeV microbursts: Evidence for dual ULF and VLF acceleration, J. Geophys. Res., 108(A8), 1329, doi:10.1029/2002JA009784, 2003.

Shumko, M., Johnson, A. T., Sample, J. G., Griffith, B. A., Turner, D. L., O'Brien, T. P., et al. (2020). Electron microburst size distribution derived with AeroCube-6. Journal of Geophysical Research: Space Physics, 125, e2019JA027651. https://doi.org/10.1029/2019JA027651