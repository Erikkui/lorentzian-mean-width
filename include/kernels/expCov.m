function C = expCov( theta, X1, X2, ~)

    sigma_f = theta(1);
    l = theta( 2);
    sigma = theta( 3);
    
    Nx1 = length( X1);
    Nx2 = length( X2);

    if Nx1 == Nx2

        C = zeros( Nx1);

        for ii = 1 : Nx1
            for jj = 1 : ii
    
                C_ij = sigma_f^2*exp( -abs( X1( ii) - X2( jj)) / l);
    
                if ii == jj
                    C_ij = C_ij + sigma^2;
                    C( ii, jj) = C_ij;
                else
                    C( ii, jj) = C_ij;
                    C( jj, ii) = C_ij;
                end
            end
        end  

    else

        C = zeros( Nx1, Nx2);

        for ii = 1 : Nx1
            for jj = 1 : Nx2
            
                C_ij = sigma_f^2*exp( -abs( X1( ii) - X2( jj)) / l);

                C( ii, jj) = C_ij;

            end  
        end

    end
end




