clc
close all
clear

includeFolders = genpath('include');
addpath( includeFolders );

%%%%%%%%%%%%%%%%% Spectrum options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bandType = 'voigt';

Nbands = 6;
Nspectrum = 1;

areaLimits = [ 1, 30];          % Peak area
gammaLimits = [ 2, 20];       % Lorentzian width, FWHM = 2*gamma
sigmaLimits = [ 2, 15];        % Gaussian width, FWHM = 2*sigma*sqrt( 2*ln(2) )
voigtMu = log(25);             % Voigt mean       
voigtSigma = 0.4;               % Voigt std
nu0Limits = [ 1625, 1675];      % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spectrumOptions.dnu = 1;
spectrumOptions.nuEdge = 35;
spectrumOptions.Nbands = Nbands;
spectrumOptions.Nsim = Nspectrum;
spectrumOptions.ALim = areaLimits;
spectrumOptions.lineWidthLimits = [ gammaLimits; sigmaLimits];
spectrumOptions.voigtParams = [ voigtMu, voigtSigma];
spectrumOptions.nu0Lim = nu0Limits;
spectrumOptions.bandType = bandType;

spectrumDataStruct = generateSpectrumData( spectrumOptions);

XData = spectrumDataStruct.nuData(:);
YData = spectrumDataStruct.noisySpectrum(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MCMC options, spectrum fit
spectrumMCMCoptions.chainLength = 20000;
spectrumMCMCoptions.covFun = @squareExpCov;

% MCMC options, FT fit
fftMCMCoptions.chainLength = 20000;
fftMCMCoptions.covFun = @squareExpCov;

MCMCoptions.spectrum = spectrumMCMCoptions;
MCMCoptions.ft = fftMCMCoptions;

outputObject = estimateMeanWidth( XData, YData, MCMCoptions);

outputObject.bandData = spectrumDataStruct.bandData;
outputObject.trueMeanGamma = spectrumDataStruct.trueMeanGamma;
%%
plotResults( outputObject)

