function [predCov, predMu] = predictFFT( theta, data, options)

    XData = data.x(:);          
    YData = data.y(:);      
    XPredict = data.XPredict(:);
    
    genCovMat = options.covFun;
    meanFun = options.meanFun;
    meanFunDerivative = options.meanFunDerivative;
    
    % Constructing Cov'(x*, x*)
    dCov_predX_predX = genCovMat( theta, XPredict, XPredict);    
    dCov_predX_predX = dCov_predX_predX - theta(3)^2*eye( size( dCov_predX_predX ) );  

    [dCov_01, dCov_10, dCov_11] = sqExpCovDerivative( theta, XPredict, XPredict);

    dCov_predX_predX = [ dCov_predX_predX, dCov_01; 
                             dCov_10,      dCov_11];    
    
    % Constructing Cov'(x*, x)
    dC_predX_X_1 = genCovMat( theta, XPredict, XData);
    [ ~, dC_predX_X_2] = sqExpCovDerivative( theta, XPredict, XData); 

    dCov_XpredX = [ dC_predX_X_1; dC_predX_X_2];             
    
    % Cov'(x, x)
    dCov_X_X = genCovMat( theta, XData, XData); 
    
    % mu'*
    meanFun1 = meanFun( theta, XData);                        
    meanFun2_1 = meanFun( theta, XPredict);
    meanFun2_2 = meanFunDerivative( theta, XPredict);
    meanFun2 = [ meanFun2_1(:); meanFun2_2(:)];              
    
    invdCov_X_X = dCov_X_X\eye( size( dCov_X_X));

    predMu = dCov_XpredX*invdCov_X_X*( YData - meanFun1) + meanFun2;
    
    % Cov'*
    predCov = dCov_predX_predX - dCov_XpredX*invdCov_X_X*dCov_XpredX';
    
    

end