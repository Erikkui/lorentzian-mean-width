function SSfun = genSSFun( theta, data, covFun, meanFun)
% Return negloglike
    
    XData = data.xdata; 
    XData = XData(:);

    YData = data.ydata; 
    YData = YData(:);
    
    N_xN_s = length( XData);

    covMat = covFun( theta, XData, XData);
    
    meanVec = meanFun( theta, XData);

    try 
    
        invCovMat = covMat\eye( size( covMat));
    
        L = chol( covMat);
        logCovDet = 2*sum( log( diag(L)));
    
        p1 = ( YData - meanVec)'*invCovMat*( YData - meanVec);
        p2 = logCovDet + N_xN_s*log( 2*pi);
    
        SSfun = p1 + p2;

    catch

        SSfun = Inf;

    end
    
end