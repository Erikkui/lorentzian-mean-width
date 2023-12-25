function [predCov, predMu] = predictSpectrum( theta, data, options)

    XData = data.x;
    YData = data.y;
    XPredict = data.XPredict;
    
    genCovMat = options.covFun;
    meanFun = options.meanFun;
     
    % Cov(x, x)
    cov_X_X = genCovMat( theta, XData, XData);
    
    % Cov(x*, x)
    cov_predX_X = genCovMat( theta, XPredict, XData);

    if length( XPredict) == length(XData)
        cov_predX_X = cov_predX_X - theta(3)^2*eye( size( cov_predX_X ) );
    end
    
    % Cov(x*, x*)
    cov_predX_predX = genCovMat( theta, XPredict, XPredict);
    cov_predX_predX = cov_predX_predX - theta(3)^2*eye( size( cov_predX_predX ) );
    
    % Predicted mean and covariance
    invCov_X_X = cov_X_X\eye( size( cov_X_X) );
    XMean = meanFun( theta, XData);
    XPredictMean = meanFun( theta, XPredict);
    
    predMu = cov_predX_X*invCov_X_X*( YData - XMean ) + XPredictMean;
    
    predCov = cov_predX_predX - cov_predX_X*invCov_X_X*cov_predX_X';
    
end

