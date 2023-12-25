function [C_01, C_10, C_11] = sqExpCovDerivative(theta, X1, X2)
%SQEXPCOVDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

sigma_f = theta(1);
l = theta( 2);

Nx1 = length( X1);
Nx2 = length( X2);

C_01 = zeros( Nx1, Nx2);
C_10 = zeros( Nx1, Nx2);
C_11 = zeros( Nx1, Nx2);
    
for ii = 1 : Nx1
    for jj = 1 : Nx2
    
        C_ij = sigma_f^2/l^2*( X1(ii) - X2(jj))*...
               exp( -0.5*( X1( ii) - X2( jj))^2 / l^2);

        C_01( ii, jj) = C_ij;

    end  
end
    
for ii = 1 : Nx1
    for jj = 1 : Nx2
    
        C_ij = sigma_f^2/l^2*( X2(jj) - X1(ii))*...
               exp( -0.5*( X1(ii) - X2(jj))^2 / l^2);

        C_10( ii, jj) = C_ij;

    end  
end   
    
for ii = 1 : Nx1
    for jj = 1 : Nx2
    
        C_ij = sigma_f^2/l^4*( l^2 - ( X1(ii) - X2(jj) )^2)*...
               exp( -0.5*( X1( ii) - X2( jj))^2 / l^2);

        C_11( ii, jj) = C_ij;

    end  
end

end

