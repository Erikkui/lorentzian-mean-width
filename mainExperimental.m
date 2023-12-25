clc
close all
clear

includeFolders = genpath('include');
addpath( includeFolders );

load data\betaCaroteneData_1240-1060.mat
% load data\betaCaroteneData_1300-1060.mat
% load data\betaCaroteneData_1360-900.mat
% load data\betaCaroteneData_1615-900.mat
% load data\betaCaroteneData_1615-1060.mat
% load data\betaCaroteneData_1615-1360.mat

XData = x;
YData = y;

figure
plot( x, y)

% MCMC options, spectrum fit
spectrumMCMCoptions.chainLength = 75000;
spectrumMCMCoptions.covFun = @squareExpCov;

% MCMC options, FT fit
fftMCMCoptions.chainLength = 75000;
fftMCMCoptions.covFun = @squareExpCov;

MCMCoptions.spectrum = spectrumMCMCoptions;
MCMCoptions.ft = fftMCMCoptions;

outputObject = estimateMeanWidth( XData, YData, MCMCoptions);

plotResults( outputObject)



