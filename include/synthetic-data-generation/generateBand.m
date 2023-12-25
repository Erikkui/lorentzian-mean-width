function [nuData, bands, trueMeanGamma] = generateBand( options)
%GENERATEBAND Generate Lorentzian, Gaussian (TODO or Voigtian) band shapes
%   type = type of band ('lorentz', 'gauss' (TODO or 'voigt')
%   Xparams = vector containing horizontal axis parameters
%   A = vector containing the band areas
%   gamma = vector containing FWHM's
%   nu0 = vector containing centers of the bands

Nbands = options.Nbands;
type = options.bandType;
ALim = options.ALim;
nu0Lim = options.nu0Lim;
lineWidthLimits = options.lineWidthLimits;
voigtParams = options.voigtParams;

deltaNu = options.dnu;
nuMin = nu0Lim(1) - options.nuEdge;
nuMax = nu0Lim(2) + options.nuEdge;

A =  randi( ALim, [Nbands,1]);       
nu0 = randi( nu0Lim, [Nbands, 1]);          

nuData = (nuMin : deltaNu : nuMax)';
dataLength = length( nuData);       

% Making the length of data even
if rem( dataLength, 2) ~= 0 
    nuData = nuData( 1 : end-1);
    dataLength = dataLength - 1;
end

% Calculating the band shapes
bands = zeros( dataLength, Nbands);
switch type
    case 'lorentz'
    
        FWHMLims = 2*lineWidthLimits( 1, :);
        FWHM = FWHMLims(1) + ( FWHMLims(2) - FWHMLims(1) )*rand( Nbands, 1);
    
        for ii = 1:Nbands
    
            nominator = A( ii)*FWHM( ii);
            denominator = pi*( FWHM(ii)^2 + 4*( nuData - nu0( ii)).^2);
            bands( :, ii) = nominator ./ denominator;
    
        end
    
        trueMeanGamma = round( mean( FWHM ), 2 );

    case 'gauss'

        sigmaLim = lineWidthLimits( 2, :);
        sigma = sigmaLim(1) + ( sigmaLim(2) - sigmaLim(1) )*rand( Nbands, 1);
        sigma = round( sigma, 2);
    
        for ii = 1:Nbands
    
            coeff = A(ii)/( sigma(ii)*sqrt( 2*pi) );
            exponent = -0.5*( nuData - nu0(ii) ).^2 ./ sigma(ii)^2;
            bands( :, ii) = coeff*exp( exponent);
        
        end
    
        trueMeanGamma = 0;

    case'voigt'   

        voigtMu = voigtParams(1) - voigtParams(2)^2/2;
        voigtSigma = voigtParams(2);

        f_V = lognrnd( voigtMu, voigtSigma, [ Nbands, 1]);  
        
        f_L_all = zeros( Nbands, 1);
        for ii = 1:Nbands
            
            f_L = f_V(ii)*rand;
            f_G = sqrt( ( f_V(ii) - f_L/2)^2 - 0.25*f_L^2 ); 
            f = ( f_G^5 + 2.69269*f_G^4*f_L + 2.42843*f_G^3*f_L^2 + ...
                 4.47163*f_G^2*f_L^3 + 0.07842*f_G*f_L^4 + f_L^5 )^(1/5);

	        %gamma = f_L/2;
	        sigma = f/( 2*sqrt( 2*log(2) ) );

            eta = 1.36603*(f_L/f) - 0.47719*(f_L/f)^2 + 0.11116*(f_L/f)^3;

            % Lorentz part
            nominator = 2*A( ii)*f;
            denominator = pi*( f^2 + 4*( nuData - nu0( ii)).^2);
            lorentz = nominator./denominator;

            % Gauss part
            coeff = A(ii)/( sigma*sqrt( 2*pi) );
            exponent = -0.5*( nuData - nu0(ii) ).^2 ./ sigma^2;
            gauss = coeff*exp( exponent); 

            % Voigt profile
            bands( :, ii) = eta*lorentz + ( 1 - eta)*gauss;

            f_L_all(ii) = f_L;
          
        end
    
        trueMeanGamma = round( mean( f_L_all ), 2 );

    otherwise

        error(['Unknown type of band shape. Band shape must be "lorentz",', ...
           ' "gauss" or "voigt".)']);

end

end
