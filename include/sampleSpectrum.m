function outputObject = sampleSpectrum( x, y, modelOptions)

    xDataInds = 1:length( y);
    dataLength = length( xDataInds);
    
    % MCMC data
    MCMC_YData = y(:);    
    MCMC_XData =  xDataInds' - 1;

    MCMCdata.xdata = MCMC_XData;
    MCMCdata.ydata = MCMC_YData;
    
    % MCMC Model
    meanFun = @( theta, X) theta(4)*ones( size( X) );
    covFun = modelOptions.covFun;
    MCMCmodel.ssfun = @( theta, data) genSSFun( theta, data, covFun, meanFun);
    MCMCmodel.N = length( MCMC_XData);
    MCMCmodel.sigma2 = 1;     
    
    % MCMC options
    MCMCoptions.nsimu = modelOptions.chainLength;
    MCMCoptions.method = 'dram';

    % MCMC parameters
    % Initial values for parameters
    sigma_S_init = std( y );
    l_init = 0.05*max( 0:size( y, 1) );   
    sigma_eps_init = sigma_S_init;
    alpha_init = mean( y, 1 );    
    MCMCparams = {
                {'\sigma_f', sigma_S_init, 0, Inf}
                {'l', l_init, 0, 2*( x(end) - x(1) )}
                {'\sigma_d', sigma_eps_init, 0, Inf}
                {'\alpha', alpha_init, 0, Inf}
                };
      
    % Running the MCMC
    [result, chain] = mcmcrun( MCMCmodel, MCMCdata, MCMCparams, MCMCoptions);
    
    %figure
    %mcmcplot( chain, [], result)
    
    % Drawing realizations of spectra
    XPredict = ( 0:length( x) - 1)';
    XPredict = [ -flip( XPredict); XPredict( 2:end)];

    theta_mean = mean( chain( end-1000:end, :), 1 );
    theta_MAP = fminsearch( MCMCmodel.ssfun, theta_mean, [], MCMCdata);

    samplingData.XPredict = XPredict;
    samplingData.x = MCMC_XData;
    samplingData.y = y;
    samplingData.Nx = dataLength;

    samplingOptions.meanFun = meanFun;
    samplingOptions.covFun = modelOptions.covFun;
    samplingOptions.theta_MAP = theta_MAP;
    samplingOptions.Nsamples = 7;
  
    sampledSpectraObj = generateSamples( samplingData, chain, samplingOptions);

    outputObject.spectra = sampledSpectraObj.sampleData;
    outputObject.theta_MAP_spectra = theta_MAP;
    outputObject.predictions = sampledSpectraObj.predictions;

end