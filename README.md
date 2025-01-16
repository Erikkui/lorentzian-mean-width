# Mean Lorentzian bandwidth estimation via two stage Gaussian Process fitting
Implementation of mean Lorentzian bandwidth estimation of a spectrum via two-stage Gaussian Process fitting, as described in [this](https://doi.org/10.1016/j.chemolab.2024.105307) article.

This method uses the fitting of two Gaussian Processes to estimate the mean Lorentzian width of a spectrum, adding uncertainty qualification for the results found in ([doi.org/10.1366/000370208784658129](https://doi.org/10.1366/000370208784658129)). The method also automates the process, eliminating the need to fit a polynomial to the mean width function described in the article. 

This is achieved by fitting a Gaussian Process (GP) to the spectrum under inspection, drawing samples of the spectrum, and fitting another GP to the Fourier Transform (FT) of the sampled spectrum. The second GP is simultaneously fitted also to the derivative of the FT, making it possible to draw samples of both the FT and its derivative. These samples can then be used to calculate the mean width and obtain a posterior distribution for it. The GP's are fitted using MCMC DRAM algorithm, which is discussed in detail [here](https://doi.org/10.1007/s11222-006-9438-0).

# Installation
Clone or download the repository and add the include folder and subfolders to your MATLAB path. To generate the synthetic data, the Statistics and Machine learning toolbox is also required.
