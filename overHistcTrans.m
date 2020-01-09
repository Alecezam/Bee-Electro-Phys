function [ histH ] = overHistcTrans( vect1, vect2, binSize, varargin )
%OVERHISTCTRANS [ histH ] = overHistcTrans( vect1, vect2, binSize, varargin )
%   Detailed explanation goes here

% Optional arguments
legLoc = 'NorthEast';
barAlpha = .5;
barWidthSetting = .85; % Scales width. A setting of 1 indicates half the bin width, I think.
%spaceBWBars = .1; %Empty space between bars between bins; 0 to 1.
overMode = 'overlay'; % Default 'overlay' is to plot bars partially 
% overlapping and transparent, but also built a mode to plot traditionally
% if you specify 'bAndw' here
g1EColor = 'r';
g2EColor = 'b';
g1MColor = 'r';
g2MColor = 'b';
statTest = 'ranksum';

pvpmod(varargin);

if strcmpi(overMode, 'bAndw') % to set up for davis lab style black and white plots
    barAlpha = 1;
    g1EColor = 'k';
    g2EColor = 'k';
    g1MColor = 'k';
    g2MColor = 'w';
    barWidthSetting = .4;
    % The offset tells the center of the group 1 bar where to go - for example,
    % offset=0.5*binSize puts the bar in the very center of the bin. Since
    % I want the control bar to have its right edge in the center of the
    % bin (i.e. at 0.5*binSize), and the default is for it to be at
    % (0.5*barWidthSetting*binSize), so I solved for
    % .5*barWidthSetting*binSize + offset = .5*binSize
    offset= 0.5 * (binSize - barWidthSetting * binSize);%0.5*binSize;
    spacer = barWidthSetting*binSize; % Moves the group 2 bar so it starts one 
    % bar-width to the right of the start of the control bar.
elseif strcmpi(overMode, 'overlay')
    offset = 0.5*barWidthSetting*binSize; %Because by default the bars are centered on the tick, and I want them between them.
    spacer = (1-barWidthSetting) * binSize;
else
    errorMess = sprintf('From overHistcTrans.m: overMode is set to [%s]. It must be the string [overlay] or [bAndw]', overMode);
    error(errorMess)
end


% Make sure the inputs are vectors and don't contain NaNs, which will break
% the mean functions.
if ~isvector(vect1)
    error('Input vector 1 is not a vector; you must enter a vector.')
end

if ~isvector(vect2)
    error('Input vector 2 is not a vector; you must enter a vector.')
end

if sum(   isnan(vect1)  )
    disp('WARNING from overhistcTrans: the first input vector contains NaN values. These will not appear in the histogram, will break the mean calculation, and will be ignored by the significance test. Use nanFix() to zero or drop them.')
end

if sum(   isnan(vect2)  )
    disp('WARNING from overhistcTrans: the second input vector contains NaN values. These will not appear in the histogram, will break the mean calculation, and will be ignored by the significance test. Use nanFix() to zero or drop them.')
end

% Generate the edges vector
xMin = min([min(vect1) min(vect2)]);
xMax = max([max(vect1) max(vect2)]);
% xLimit(1) = round2(xMin-binSize,binSize);
% xLimit(2) = round2(xMax+binSize,binSize);
xLimit(1) = (xMin-binSize);
xLimit(2) = (xMax+binSize);
inEdges = xLimit(1):binSize:xLimit(2);


% Find means of each input vector
mu1 = mean(vect1);
mu2 = mean(vect2);

% Bin vectors into histograms
hist1 = histc(vect1, inEdges);
n1 = sum(hist1);
hist1freq = hist1; %%Edited by alec to make it in counts instead of proportion
hist2 = histc(vect2, inEdges);
n2 = sum(hist2);
hist2freq = hist2;

% Make the figure
histH = figure;



bars1H = bar(inEdges + offset, hist1freq, g1MColor, 'BarWidth', barWidthSetting);
%bars1H = bar(inEdges -barWidthSetting*binSize, hist1freq, 'k', 'BarWidth', barWidthSetting);
hold on
bars2H = bar(inEdges + offset + spacer, hist2freq, g2MColor, 'BarWidth', barWidthSetting, 'LineWidth', 1);
%bars2H = bar(inEdges, hist2freq, 'r', 'BarWidth', barWidthSetting);

% Get max value of histograms for plotting asterisks for means
yMax = max([max(hist1freq) max(hist2freq)]);
mu1LabH = scatter(mu1, yMax+.025, 100, 'MarkerEdgeColor', g1EColor, 'MarkerFaceColor', g1MColor);
mu2LabH = scatter(mu2, yMax+.025, 100, 'MarkerEdgeColor', g2EColor, 'MarkerFaceColor', g2MColor);
set(gca, 'FontSize', 18)
legend([bars1H bars2H], 'Linalool', 'Ethanol', 'Location', legLoc)
ylabel('Number of Spikes');
titleStr = sprintf('%d Linalool; %d Ethanol' ,n1, n2);
title(titleStr);

set(gca, 'XLim', [xLimit(1)-binSize xLimit(2)+binSize]);
set(gca, 'xtick', inEdges);
set(gca, 'YLim', [0 yMax+.1])

condense = ceil(length(inEdges)/7); % factor by which labels must be condensed
xTickLabs = cell(1,length(inEdges));

% Set tick labels so that 8 are shown, including first and last.
for i = 1:length(inEdges)
    if ~mod( (i-1) , condense)
     xTickLabs{i} = num2str(inEdges(i));
    end
end
xTickLabs{end} = num2str(inEdges(end));

set(gca, 'XTickLabel', xTickLabs)% % XTickLabel

% Test for signicicance
if strcmpi(statTest, 'ranksum')
    [p, h, stats] = ranksum(vect1, vect2, 'alpha', 0.05, 'tail', 'both');
elseif strcmpi(statTest, 'ttest2')
    [h p] = ttest2(vect1, vect2, .05, 'both');
else
    error('From overHistcTrans.m: statTest (an optional argument) is not recognized. Must be the string ranksum or ttest2  .')
end



%%
%{
%Trying to find if any of the bins are significantly sized
nBins = size(hist1)
nBins = nBins{1};
popMean = n1 /nBins;

popSTD = 0;
for i = 1: nBins
    popSTD = popSTD + (hist1(i) - popMean) * (hist1(i) - popMean);
end
popSTD = popSTD /n1;
popSTD = sqrt(popSTD);

zScores = {};
for i = 1:nBins
    z = (hist1(i) - popMean)/(popSTD/sqrt(n1));
    zScores = [zScores; {z, 'linalool'}];
end


nBins = size(hist2)
nBins = nBins{1};
popMean = n2/ nBins;

popSTD = 0;
for i = 1: nBins
    popSTD = popSTD + (hist2(i) - popMean) * (hist2(i) - popMean);
end
popSTD = popSTD /n2;
popSTD = sqrt(popSTD);

zScores = {};
for i = 1:nBins
    z = (hist2(i) - popMean)/(popSTD/sqrt(n2));
    zScores = [zScores; {z, 'ethanol'}];
end

newFile = strcat('C:\Users\Alec\Desktop\hokanson lab bee elecphys\ephys test data\Outputs', '\zScores');

print(newFile, '-dpng');

%}
%%

% Draw line between groups
lineLevel = yMax + .035;
lineH = line([mu1 mu2], [lineLevel lineLevel]);
set(lineH, 'Color', [0 0 0], 'LineWidth', 2);
% % Make a label for the line
% % Determine significance level
if p < .001
    sig = '***';
elseif p<.01
    sig = '**';
elseif p<.05
    sig = '*';
else
    sig = 'n.s.';
end

% 
% % Create text string
% lineNote = sprintf('%s, p=%1.3f', sig, p);
% % Make the label
% hSig = text((mu1+mu2)/2, lineLevel+.015, lineNote, 'HorizontalAlignment', 'center');
% set(hSig, 'FontSize', 14);

alpha(barAlpha)

end

