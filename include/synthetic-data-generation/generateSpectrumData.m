function data = generateSpectrumData( options)
% GENERATESPECTRUMDATA Generate multiple spectra consisting of defined
% number of bands 
%   Nbands = number of bands of which the spectra is built of
%   Nsim = number of spectra to be generated
%   Options = struct containing the limits of x-axis and range of band 
%             areas (A), FWHM's (gamma) and band centers (nu0)  

Nsim = options.Nsim;
Fs = 1/options.dnu;                 % Sample rate

%%%%%%%%%%% Generating the bands and the spectrum

% Lorentzian band(s)
[nuData, bands, trueMeanGamma] = generateBand( options);

figure
bandPlot = plot( nuData, bands, 'b-');
title( 'Spectrum without noise and bands')

% Calculating the spectrum
spectrum = sum( bands, 2);

hold on
grid on
spectrumPlot = plot( nuData, spectrum, 'k-', 'LineWidth', 1);
legend( [ bandPlot( end), spectrumPlot], 'Bands', 'Spectrum')


%%%%%%%%%%%% Calculating the mean width           
dataLength = length( spectrum);
noiseSigma = 0.05*max( spectrum);

noisySpectrumData = zeros( dataLength, Nsim);
fftSpectrumData = zeros( dataLength/2 + 1, Nsim);
fftDerivativeData = zeros( dataLength/2 - 1, Nsim);
meanGammaData = fftDerivativeData;

xData = ( Fs*( 0:(dataLength/2) )/dataLength )';
deltaX = abs( xData( 2) - xData(1)); 

for ii = 1:Nsim

    noisySpectrum_temp = spectrum + noiseSigma*randn( dataLength, 1);
    noisySpectrumData( :, ii) = noisySpectrum_temp;
    
    % FFT 
    [ fftSpectrum_temp, ~] = onesidefft( noisySpectrum_temp, Fs);
    fftSpectrumData( :, ii) = fftSpectrum_temp;
    
    % FFT Derivatives, central difference
    fftDerivative_temp = ( fftSpectrum_temp( 3:end) - fftSpectrum_temp( 1:end-2) );
    fftDerivative_temp = fftDerivative_temp/( 2*deltaX);
    fftDerivativeData( :, ii) = fftDerivative_temp;
    
    %Lorentzian width function
    meanWidth_temp = -fftDerivative_temp ./ ( pi*fftSpectrum_temp( 2:end-1) );
    meanGammaData( :, ii) = meanWidth_temp;
end

data.nuData = nuData;
data.xData = xData;
data.spectrum = spectrum;
data.noisySpectrum = noisySpectrumData;
data.bandData = bands;
data.fftSpectrum = fftSpectrumData;
data.fftDerivative = fftDerivativeData;
data.meanGamma = meanGammaData;
data.noiseSigma = noiseSigma;
data.trueMeanGamma = trueMeanGamma;

end