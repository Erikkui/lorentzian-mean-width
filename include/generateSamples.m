function outputObject = generateSamples( data, chain, options)
    
    Nx = data.Nx;
    NxPred = 2*Nx - 1;
	Nsamples = options.Nsamples;
    chainLength = size( chain, 2);

    burnedChain = chain( round( chainLength/2 ):end, :);
    thetaInds = randi( length( burnedChain), Nsamples, 1);
    sampleThetas = burnedChain( thetaInds, :);
  
    if isfield( options, 'meanFunDerivative')
        
        Nrand = options.Nrand;
        indsFFT = Nx : 2*Nx - 1;
        indsFFTDer = 3*Nx - 1 : 4*Nx - 2;
        
        sampleData = zeros( Nrand*2*NxPred, Nsamples);          
        for ii = 1:Nsamples
    
            theta = sampleThetas( ii, :);
            [predCov, predMu] = predictFFT( theta, data, options );
            
            predCovFFTDiag = sqrt( diag( predCov ));
            fftMin = min( abs( predCovFFTDiag( indsFFT) ) );
            dfftMin = min( abs( predCovFFTDiag( indsFFTDer ) ) );

            nugget = 1e-4*fftMin*ones( 2*Nx-1, 1);
            dNugget = 1e-4*dfftMin*ones( 2*Nx-1, 1);
            totalNugget = diag( [ nugget; dNugget]); 

            samplingPredCov = predCov + totalNugget;

            sampleData( :, ii) = cholSample( predMu, samplingPredCov, Nrand);

        end     

        sampleData = reshape( sampleData, 2*NxPred, Nrand*Nsamples);

        sampledFFT = sampleData( indsFFT, :);
        sampledFFTDerivative = sampleData( indsFFTDer, :);
        [predCov, predMu] = predictFFT( options.theta_MAP, data, options );

        outputObject.sampledFFT = sampledFFT;
        outputObject.sampledFFTDerivative = sampledFFTDerivative;
        outputObject.predictions = struct( 'predMu', predMu, 'predCov', predCov);

    else

        sampledData = zeros( Nx, Nsamples);
        for ii = 1:Nsamples
    
            theta = sampleThetas( ii, :);
            [predCov, predMu] = predictSpectrum( theta, data, options );
            
            samplingPredCov = predCov + 1e-4*eye( size( predCov ) );

            sampleData = cholSample( predMu, samplingPredCov, 1);
            sampleData = sampleData + theta(3)*randn( size( sampleData) );
            
            sampledData( :, ii) = sampleData( Nx:end);
            
        end

        [predCov, predMu] = predictSpectrum( options.theta_MAP, data, options );
        predMu = predMu( Nx:end);
        predCov = predCov( Nx:end);

        outputObject.sampleData = sampledData;
        outputObject.predictions = struct( 'predMu', predMu, 'predCov', predCov);
    end
end