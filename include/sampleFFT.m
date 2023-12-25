function outputObject = sampleFFT( x, y, modelOptions)

    Nspect = size( y, 2);

    fftXDataInds = ( 1:30 )';

    % MCMC data
    MCMC_XData =  fftXDataInds - 1;
    MCMC_XData = repmat( MCMC_XData, Nspect, 1);

    MCMC_YData = y( fftXDataInds, :);
    MCMC_YData = MCMC_YData(:);
    
    MCMCdata.xdata = MCMC_XData;
    MCMCdata.ydata = MCMC_YData;

    % MCMC Model
    meanFun = @( theta, X) theta(4)*exp( theta(5)*X );
    covFun = modelOptions.covFun;
    MCMCmodel.ssfun = @( theta, data) genSSFun( theta, data, covFun, meanFun);
    MCMCmodel.N = length( MCMC_XData);
    MCMCmodel.sigma2 = 1;

    % MCMC options
    MCMCoptions.nsimu = modelOptions.chainLength;
    MCMCoptions.method = 'dram';

    % MCMC parameters
    % Initial values for parameters
    sigma_f_init = std( y(:) ); 
    lambda_init = 0.1*max( fftXDataInds );   
    sigmaD_init = sigma_f_init;
    alpha_init = 0.2*max( y, [], 'all');
    beta_init = 1;
    
    MCMCparams = {
                {'\sigma_f', sigma_f_init, 0, Inf}
                {'lambda', lambda_init, 0, 3*max( fftXDataInds)}
                {'\sigma', sigmaD_init, 0, Inf}
                {'\beta_1', alpha_init, 0, 10*max( y, [], 'all')}
                {'\beta_2', beta_init, -Inf, Inf}
              };
        
    % Running the MCMC
    [results, chain] = mcmcrun( MCMCmodel, MCMCdata, MCMCparams, MCMCoptions);
    
    %figure
    %mcmcplot( chain, [], results)

    % Draw realizations of fft and its derivative
    XPredict = fftXDataInds - 1;
    XPredict = [ -flip( XPredict); XPredict( 2:end)];

    theta_mean = mean( chain( end-1000:end, :), 1 );
    theta_MAP = fminsearch( MCMCmodel.ssfun, theta_mean, [], MCMCdata);

    samplingData.XPredict = XPredict;
    samplingData.x = MCMC_XData;
    samplingData.y = MCMC_YData;
    samplingData.Nx = length( fftXDataInds);
    
    samplingOptions.meanFun = meanFun;
    samplingOptions.meanFunDerivative = @( theta, X) theta(4)*theta(5)*exp( theta(5)*X );
    samplingOptions.covFun = modelOptions.covFun;
    samplingOptions.theta_MAP = theta_MAP;
    samplingOptions.Nsamples = 700;
    samplingOptions.Nrand = 5;
  
    sampledFFTObj = generateSamples( samplingData, chain, samplingOptions);
    dXi = abs( x(2) - x(1) );
    
    outputObject.sampledFFT = sampledFFTObj.sampledFFT;
    outputObject.sampledFFTDerivative = sampledFFTObj.sampledFFTDerivative / dXi;
    outputObject.theta_MAP = theta_MAP;
    outputObject.predictions = sampledFFTObj.predictions;
end