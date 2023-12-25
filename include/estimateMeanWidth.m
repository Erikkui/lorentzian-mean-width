function outputObject = estimateMeanWidth( x, y, MCMCoptions)

    dx = abs( x(2) - x(1) );
    Fs = 1/dx;
    SZ = size( y, 1);

    if rem( SZ, 2) ~= 0 
        x = x( 1 : end-1);
        y = y( 1 : end-1);
        SZ = SZ - 1;
    end

    spectrumMCMCoptions = MCMCoptions.spectrum;
    fftMCMCoptions = MCMCoptions.ft;

    spectrumSampleDataObj = sampleSpectrum( x, y, spectrumMCMCoptions);

    %%%
    sampledSpectra = spectrumSampleDataObj.spectra;
    Nspect = size( sampledSpectra, 2);
    fftX = ( Fs*( 0:(SZ/2) )/SZ )';
        
    fftY = zeros( SZ/2+1, Nspect);
    for ii = 1:Nspect
        fftY( :, ii) = onesidefft( sampledSpectra(:, ii), Fs); 
    end

    fftSampleDataObj = sampleFFT( fftX, fftY, fftMCMCoptions);

    fftSamples = fftSampleDataObj.sampledFFT;
    fftDerivativeSamples = fftSampleDataObj.sampledFFTDerivative;

    meanGamma = -fftDerivativeSamples ./ ( pi*fftSamples);

    %%% 
    
    outputObject.nu = x;
    outputObject.spectrum = y;
    outputObject.spectrumSamples = spectrumSampleDataObj;
    outputObject.fftX = fftX;
    outputObject.fftY = fftY;
    outputObject.FTSamples = fftSampleDataObj;
    outputObject.meanGamma = meanGamma;

end