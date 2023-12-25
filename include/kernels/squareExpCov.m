function C = squareExpCov( theta, X1, X2, ~)

    sigma_f = theta(1);
    l = theta( 2);
    sigma = theta( 3);
    
    Nx1 = length( X1);
    Nx2 = length( X2);

    if Nx1 == Nx2

        C = zeros( Nx1);

        for ii = 1 : Nx1
            for jj = 1 : Nx2
    
                C_ij = sigma_f^2*exp( -1/2*( X1( ii) - X2( jj))^2 / l^2);
                
                if ii == jj
                    C_ij = C_ij + sigma^2;
                    C( ii, jj) = C_ij;
                else
                    C( ii, jj) = C_ij;
                    C( jj, ii) = C_ij;
                end
            end
        end  
        % Adding nugget for numerical stability
       % C = C + eye( size( C))*1e-8;

    else

        C = zeros( Nx1, Nx2);

        for ii = 1 : Nx1
            for jj = 1 : Nx2
            
                C_ij = sigma_f^2*exp( -0.5*( X1( ii) - X2( jj))^2 / l^2);

                C( ii, jj) = C_ij;

            end  
        end

    end
end




