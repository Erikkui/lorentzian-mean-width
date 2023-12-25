function sample = cholSample( mu, cov, Nrand)  
      
    try
        L = chol( cov);
    
        normalSamples = randn( length( mu ), Nrand);
        
        sample = mu + L'*normalSamples;
    
        sample = sample(:);
    catch
        sample = NaN*ones( length(mu), Nrand );
        sample = sample(:);
    end
end