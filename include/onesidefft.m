function [YData, XData] = onesidefft(YData, Fs)
%ONESIDEFFT Calculates one-sided FFT for given (1-D) data vector
%   

SZ = size( YData, 1);
FFT_YData = fft( YData, [], 1);
P2 = abs( FFT_YData./(2*Fs) );        % Kysy teemulta miksi skaalaus 2*Fs
P1 = P2( 1:SZ/2 + 1, :);
P1( 1:end-1, :) = 2*P1( 1:end-1, :);

YData = P1;
XData = (Fs*( 0:(SZ/2))/SZ)';

end
