function plotResults( dataObj)

nuData = dataObj.nu;
spectrumData = dataObj.spectrum;
sampledSpectra = dataObj.spectrumSamples.spectra;

fftInds = 1:30;
xiData = dataObj.fftX( fftInds);
fftData = dataObj.fftY( fftInds, :);

fftSampleData = dataObj.FTSamples.sampledFFT( fftInds, :);
fftDerivativeSampleData = dataObj.FTSamples.sampledFFTDerivative( fftInds, :);

meanGamma = dataObj.meanGamma( fftInds, :);

meanGammas_0 = ( meanGamma( 1, :) );
meanGammas_0 = meanGammas_0( meanGammas_0 > 0);

%%%% Plotting spectrum and bands if exist
spectFig = figure();
hold on
axis tight

if isfield( dataObj, 'bandData')
    bandData = dataObj.bandData;     
    bandPlot = plot( nuData, bandData, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
end 

spectPlot = plot( nuData, spectrumData);
spectPlot.LineStyle = "-";
spectPlot.LineWidth = 2;
spectPlot.Color = "#0072BD";
spectPlot.DisplayName = "\boldmath$S$";

h = xlabel( '$\nu$');
h.Interpreter = 'latex';

lg = legend( spectPlot);
lg.Interpreter = 'latex';

hold off


%%%% Plotting spectrum samples
spectSampleFig = figure();
hold on
axis tight

spectSamplePlot = plot( nuData, sampledSpectra, 'Color', '#0072BD', ...
                        'LineWidth', 1.25, 'HandleVisibility', 'off');
spectSamplePlot(1).DisplayName = "{\boldmath$\widetilde{S}$}";

h = xlabel( '$\nu$');
h.Interpreter = 'latex';

lg = legend( spectSamplePlot(1) );
lg.Interpreter = 'latex';

hold off

% Setting reverse x-axis for spectrum plots
fh = findall(0, 'Type', 'Axes');
set( fh, 'xdir', 'reverse')


%%% Plotting fft of sampled spectra
fftFig = figure();
hold on
axis tight

fftDataPlot = plot( xiData, fftData(:, 1:7) , 'Color', '#0072BD', ...
                    'LineWidth', 1.5, 'HandleVisibility', 'off');
fftDataPlot(1).DisplayName = ' \boldmath$ Z$';

h = xlabel( '$\xi$');
h.Interpreter = 'latex';

lg = legend( fftDataPlot(1) );
lg.Interpreter = 'latex';
lg.Location = 'east';

hold off


% %%%% Plotting sampled fft and its derivative
% left_color = [0 0 0];
% right_color = [0 0.4470 0.7410];
% 
% fftSampleFig = figure();
% 
% set( fftSampleFig, 'defaultAxesColorOrder', [ left_color; right_color]);
% 
% hold on
% axis tight
% yyaxis right
% 
% for ii = 1:7
%     FFTDerivativeSamplePlot = plot( xiData, fftDerivativeSampleData( :, ii), ...
%                              'Color', [0 0.4470 0.7410], 'LineWidth', 1.5,  ... 
%                              'HandleVisibility', 'off');
% end
% FFTDerivativeSamplePlot.DisplayName = "{\boldmath$\widetilde{Z'}$}";
% 
% axis tight
% yyaxis left
% 
% for ii = 1:7
%     FFTsamplePlot = plot( xiData, fftSampleData( :, ii),  'Color', [0 0 0], ...
%                        'LineWidth', 1.5, 'HandleVisibility', 'off');
% end
% FFTsamplePlot.DisplayName = "{\boldmath$\widetilde{Z}$}";
% 
% hx = xlabel( '$\xi$');
% hx.Interpreter = 'latex';
% 
% lg = legend( [ FFTsamplePlot, FFTDerivativeSamplePlot ] );
% lg.Interpreter = 'latex';
% lg.Location = "east";
% 
% hold off


%%%% Plotting mean gamma
meanGammaFig = figure();
hold on
axis tight

meanGammaPlot = plot( xiData, meanGamma, 'LineWidth', 1, 'Color', '#0072BD', ...
                        'HandleVisibility', 'off');
meanGammaPlot(1).DisplayName = "$\widetilde{\gamma}(\xi)$";
lx = xlabel('$\xi$');
lx.Interpreter = 'latex';

lg = legend( meanGammaPlot(1) );
lg.Interpreter = 'latex';
lg.Location = 'northwest';

hold off


%%%% Histogram of mean gamma at xi = 0
meanGammaHist = figure();
hold on
axis tight

h = histogram( meanGammas_0);
h.Normalization = "probability";
h.DisplayName = '$p(\overline{\gamma}\mid${\boldmath$Z$}$)$';

if isfield( dataObj, 'trueMeanGamma')
    trueMeanGamma = dataObj.trueMeanGamma;
    xl = xline( trueMeanGamma, 'LineWidth', 8, 'Color', 'black');
    xl.HandleVisibility = 'off';
end

lx = xlabel('$\overline{\gamma}$');
lx.Interpreter = 'latex';

lg = legend;
lg.Interpreter = 'latex';

hold off


%%% Setting font sizes
fh = findobj('-property', 'FontName');
set( fh, 'FontSize', 26.5)
fh = findall(0, 'Type', 'Axes');
set( fh, 'TitleFontSizeMultiplier', 1.5)

end