clear
clc
tic

headerFile = uigetdir('','Select header file');
fprintf("Header: %s\n", headerFile);
inFile = dir (headerFile);
dateFiles = {};
numDateFiles = 0;
for i = 3:size(inFile)-2 %Starting at 3 bc there are 2 files of nonsense
    dateFiles = [dateFiles; {inFile(i).name}]; %add a new row of data
    numDateFiles = numDateFiles + 1;
end

hiveSessionFiles = {}; %This will store paths to all of the session files for each bee
for i = 1:numDateFiles
    holderStr = strcat(headerFile,'\',dateFiles(i));
    holderDir = dir(holderStr{1});
    for j=1:size(holderDir)
       if(~strcmp(holderDir(j).name,'.') & ~strcmp(holderDir(j).name,'..') & ~strcmp(holderDir(j).name,'.DS_Store'))
           holderPath = strcat(holderDir(j).folder,'\',holderDir(j).name);
           hiveSessionFiles = [hiveSessionFiles; {holderPath, holderDir(j).name}];
       end
    end
end

experimentFiles = dir(strcat(headerFile,'\Experiments'));

hiveFiles = {}; %This will store the path to all of the experiment hive files

for i = 1:size(experimentFiles)
    if(~strcmp(experimentFiles(i).name,'.') & ~strcmp(experimentFiles(i).name,'..') & ~strcmp(experimentFiles(i).name,'.DS_Store'))
        holderPath = strcat(experimentFiles(i).folder,'\',experimentFiles(i).name);
        hiveFiles = [hiveFiles; {holderPath, experimentFiles(i).name}];
    end
end

outputFile = strcat(headerFile,'\Outputs');
for k = 1:size(hiveSessionFiles)
   %Starting with data in 062118 folder
%b03 = bee 3, meaning from the third hive
%s01 = session 1; the first bee from the hive mentioned before
% The experiments folder has the odor presentations etc (.mat files).
% "bee_001" contains all the sessions for hive 1

% I believe Priya's data is saved in "One File Per Channel" format

% Issues: lots of 60 Hz noise. How big are signals supposed to be? The 60hz
% noise is significant but seems small because some other signals are huge.
% Are they real?

% Need to add the Google Drive Electrophys Core/MATLAB folder and
% subfolders to directory before running this


%% Create data vectors from Priya's data


%HeaderDataDir = '/Users/hokansok/Google Drive/OSU/Electrophys Core/MATLAB/ephys test data/062518/b06s01_180625_173911';
%ExperimentInfoDir = '/Users/hokansok/Google Drive/OSU/Electrophys Core/MATLAB/ephys test data/Experiments/bee_006';
%ExperimentInfoMatfile = 'm006s01.mat';

%%

HeaderDataDir = hiveSessionFiles(k);
HeaderDataDir = HeaderDataDir{1};
for i=1:size(hiveFiles)
    currentExp = strcat('bee_0',extractBetween(hiveSessionFiles{k,2},2,3));
    currentExp = currentExp{1};
    if(hiveFiles{i,2} == currentExp)
        ExperimentInfoDir = hiveFiles(i);
        ExperimentInfoDir = ExperimentInfoDir{1};
    end
end

ExperimentInfoMatfile = strcat('m0',extractBetween(hiveSessionFiles{k,2},2,6),'.mat'); %'m006s01.mat';
ExperimentInfoMatfile = ExperimentInfoMatfile{1};

display(ExperimentInfoMatfile);
display(hiveSessionFiles{k,2});

startOfDate = (size(headerFile)+ 1);
startOfDate = startOfDate(2);
endOfDate = (size(headerFile)+8);
endOfDate = endOfDate(2);

beeOutputFile = strcat(outputFile, extractBetween(hiveSessionFiles{k,1},startOfDate,endOfDate));
beeOutputFile = strcat(beeOutputFile,ExperimentInfoDir(end-2:end),extractBetween(hiveSessionFiles{k,2},4,6)); 
beeOutputFile = beeOutputFile{1};
status = mkdir(beeOutputFile);

%fid = fopen(HeaderDataDir, 'r');

%%
cd(HeaderDataDir);
LowpassFiltCutoff = 1000; % Default: 1 kHz

disp('Select the corresponding RHD (header file): ephys test data/062518/b06s01_180625_173911/info.rhd')
% Read in the header file
read_Intan_RHD2000_file; %select the corresponding RHD (header file) 
%Above file slightly edited by Alec to auto read info.rhd file based on
%current path
% Pull sampling rate from the parameters structure imported by header file
sr = frequency_parameters.amplifier_sample_rate; %Hz

% Extract time vector (size: numSamples x 1) in seconds, from the time.dat file
fileinfoTimedat = dir('time.dat');
num_samples_time = fileinfoTimedat.bytes/4; % int32 = 4 bytes ; this value is hardcoded in the instruction pdf
fid = fopen('time.dat', 'r');
timeVector = fread(fid, num_samples_time, 'int32'); % this value is hardcoded in the instruction pdf
fclose(fid);
timeVector = timeVector / sr; % sample rate from header file

% Extract the data from an amplifier channel .dat file
% Data is saved in v, a vector of size numSamples x 1, in microvolts
fileinfoAmp017 = dir('amp-B-017.dat');
num_samples_amp17 = fileinfoAmp017.bytes/2; % int16 = 2 bytes ; this value is hardcoded in the instruction pdf
fid = fopen('amp-B-017.dat', 'r');
v = fread(fid, num_samples_amp17, 'int16'); % this value is hardcoded in the instruction pdf
fclose(fid);
v = v * 0.195; % convert to microvolts
toc

%% Override data vectors with sample data for testing

UseTestData = 0;
if UseTestData
    disp('Warning from priyaslush: data vectors are being overridden with sample data')
    load('/Users/hokansok/Google Drive/OSU/Electrophys Core/MATLAB/spikesortbootcamp/SpikeSortData.mat')
    timeVector = Time;
    v = Vtotal;
end

%% Analyze data vectors

% Constants for figures
FWidth = 1400; % Width of sample traces
FHeight = 200; % Height of sample traces

plotUnfilt = 0;
if plotUnfilt
    fig1Handle = figure;
    plot(timeVector, v)
    
    % Format
    set(fig1Handle,'PaperPositionMode','auto')
    set(fig1Handle, 'Position', [0 600 FWidth FHeight])
end

% Set up Butterworth Lowpass filter
%What's the lowest order usable for the filter? https://www.mathworks.com/help/signal/ref/buttord.html
passbandEdge = 1000 / (sr / 2); %output is in pi rad/samples
stopbandEdge = 2000 / (sr / 2);
RippleAllowed = 3; % dB
AttenuationInStop = 60; % dB
[n,Wn] = buttord(passbandEdge,stopbandEdge,RippleAllowed,AttenuationInStop);
% design filter
[b,a] = butter(n,Wn);
% filter
vFilt = filtfilt(b, a, v);

% Create summary figure
sumFigHandle = figure;
plot(timeVector, v, 'Color', [.5 .5 .5])
hold on

% Format
set(sumFigHandle,'PaperPositionMode','auto')
set(sumFigHandle, 'Position', [0 600 FWidth FHeight])
title(strcat(hiveSessionFiles{k,2}(1:13),' - amp-B-017.dat'))
xlabel('Time (s)')
ylabel('Response (uV)')
YMax = 5000;
YMin = -5000;
ylim([YMin YMax]);
if max(abs(vFilt)) > YMax
    disp('WARNING from priyaslush: the summary figure y limits are set so that some of the trace is not being shown.')
end


% Divide into epochs
cd(ExperimentInfoDir);
load(ExperimentInfoMatfile)

% Construct odor "on" matrices [onset(ms) offset(ms) concentration(1-100)]
numParameters = size(tr,2);
i=1;
jLin=1;
jEth=1;
for i = 1:numParameters
    if strcmp(tr(i).behavPar, 'OdorPresentation')
        if strcmpi(tr(i).odorName, 'linalool')
            %Produce linalool matrix
            linaloolMat(jLin,:) = [tr(i).t0    tr(i).t0+tr(i).odorDur   tr(i).odorConc];
            jLin = jLin+1;
        elseif strcmpi(tr(i).odorName, 'ethanol')
            %Produce ethanol matrix
            ethanolMat(jEth,:) = [tr(i).t0    tr(i).t0+tr(i).odorDur   tr(i).odorConc];
            jEth = jEth+1;
        else
            i
            error('From priyaslush: when looking through the tr structure (of experimental info), I saw an odor I have never seen before. I should edit the loop to define whath to do with that odor. Index was printed just above.')
        end
    else
        i
        error('From priyaslush: when looking through the tr structure (of experimental info), I saw an experimental parameter I have never seen before. I should edit the loop to define whath to do with that parameter. Index was printed just above.')
    end
end

% Set colors for the odor concentration plots and for the optional shading
% on the summary plot
linColor = [1 0 0];
ethColor = [0 0 1];

% Generate odor concentration traces, and optionally add shading to summary
% plot


if not(UseTestData) %This code doesn't work if you use test data instead of Priya's data
    
    linaloolConcTrace = makeOdorTrace(timeVector, sr, linaloolMat, 'shadeColor', linColor);
    ethanolConcTrace = makeOdorTrace(timeVector, sr, ethanolMat, 'shadeColor', ethColor);
    ylim([YMin YMax]);
    
    fig3Handle = figure;
    linaloolh = plot(timeVector, linaloolConcTrace);
    linaloolh.Color = linColor;
    
    hold on
    ethanolh = plot(timeVector, ethanolConcTrace);
    ethanolh.Color = ethColor;
    
    
    % Format
    ylim([YMin YMax]);
    set(fig3Handle,'PaperPositionMode','auto')
    set(fig3Handle, 'Position', [0 300 FWidth FHeight])
    title(strcat('Odors bee-006-',ExperimentInfoMatfile))
    xlabel('Time (s)')
    ylabel('Odor Concentration')
end
toc
% %}

% Identify a baseline period and get stdev during it
disp('BaselinePeriod window is currently hard-coded')
StdevThreshold = 300; % What percent of the stdev should be used as the threshold?
BaselinePeriod = [475 495]; %manually chosen not-very-noisy period
% Highlight the baseline period as orange on the summary plot
baseStart = find(timeVector == BaselinePeriod(1)); %Sample of start of baseline
baseEnd = find(timeVector == BaselinePeriod(2)); %Sample of end of baseline

% Plot the part of the trace used for baseline
figure(sumFigHandle); % Make the summary figure the active figure
BaselineStdevColor = [1 .7 0];
plot(timeVector(baseStart:baseEnd), vFilt(baseStart:baseEnd), 'Color', BaselineStdevColor)
hold on
% Calculate and plot stdev thresholds
if abs(mean(vFilt(baseStart:baseEnd)))>1
   % error('Error from priyaslush: the mean value of the baseline period is more than 1 uV away from 0. Check to be sure this is ok, or consider loosening the threshold')
end

if abs(mean(vFilt(baseStart:baseEnd)) - mean(vFilt)) >1
   % error('Error from priyaslush: the mean value of the baseline period is more than 1 uV away from the mean value of the entire trace. Check to be sure this is ok, or consider loosening the threshold')
end
stDev = std(vFilt(baseStart:baseEnd));
spikeThresh = stDev * StdevThreshold/100;
stdLinesH = line( [timeVector(1) timeVector(1) ; timeVector(end) timeVector(end)] , [spikeThresh -spikeThresh ; spikeThresh -spikeThresh], 'Color', BaselineStdevColor);
hold on;
%

%%
%for finding the spikes and getting data on them
%collects start, end, peak time/value, scent type, power, and time since
%scent
%{
spikeInfo = {};

start_time = 1;
end_time = size(timeVector);
start_of_spike = 1;
end_of_spike = 1;

for i = start_time:end_time
    if(vFilt(i) > spikeThresh && start_of_spike == 1)
        start_of_spike = timeVector(i);
    elseif(vFilt(i) < spikeThresh && start_of_spike ~= 1)
        end_of_spike = timeVector(i);
            timeOfSpike = findAboveSpike(start_of_spike, end_of_spike, vFilt);
        if(timeOfSpike == 0)
            fprintf("No spike found between %5.5f and %5.5f\n", start_of_spike, end_of_spike);
            %above line should never print, there should always be a spike found
        else
            lastScent = findLastScent(timeOfSpike, linaloolConcTrace, ethanolConcTrace);
            if(strcmp(lastScent(1,1),""))
                fprintf("ERROR NO LAST SCENT FOUND");
            end
            spikeInfo = [spikeInfo; {start_of_spike, end_of_spike,timeOfSpike/sr,vFilt(timeOfSpike), lastScent{1,1}, lastScent{1,2}, lastScent{1,3}/sr}]; 
            %^^add a new row of spike data, contains time when spike passes
            %cutoff line at the start and end, the time of the local peak,
            %and the value at that time

        end
        start_of_spike = 1;
        end_of_spike = 1;
    end
end

%marks the spikes on the plot of the voltages
%red for linalool, blue for ethanool
numSpikes = size(spikeInfo);
numSpikes = numSpikes(1,1);
for i = 1:numSpikes
   % markSize = spikeInfo{i,7}/(31/6);
   % marSize = fix(markSize);
    if(strcmp(spikeInfo{i,5},"linalool"))
        plot(spikeInfo{i,3}, spikeInfo{i,4}, 'r*','MarkerSize',3);
    elseif(strcmp(spikeInfo{i,5},"ethanol"))
        plot(spikeInfo{i,3}, spikeInfo{i,4}, 'b*','MarkerSize',3);
    end
    hold on;
end

scentPower = cell2mat(spikeInfo(:,6));

linaloolScentPowers = {};
ethanolScentPowers = {};

for i = 1:numSpikes
   if(strcmp(spikeInfo{i,5},"linalool"))
        linaloolScentPowers = [linaloolScentPowers; {spikeInfo{i,6}}];
   elseif(strcmp(spikeInfo{i,5},"ethanol"))
        ethanolScentPowers = [ethanolScentPowers ; {spikeInfo{i,6}}];
   end
end

linaloolScentPowers = cell2mat(linaloolScentPowers);
ethanolScentPowers = cell2mat(ethanolScentPowers);

figure;
histogram(linaloolScentPowers, 'FaceColor', 'r', 'BinMethod', 'integers');
title(strcat(hiveSessionFiles{k,2}(1:13),' - Odor Histogram'));
ylabel("Counts");
xlabel("Odor Strength (?)");
hold on;
histogram(ethanolScentPowers, 'FaceColor', 'b');
%}


%%
%This script filters out the 60hz noise, and then finds spikes based on the
%filtered data, then plots the non-filtered data with the found spikes
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',sr);
y = filtfilt(d,v);
 
% figure;
% plot(timeVector, y);
% hold on;
% stdLinesH = line( [timeVector(1) timeVector(1) ; timeVector(end) timeVector(end)] , [spikeThresh -spikeThresh ; spikeThresh -spikeThresh], 'Color', BaselineStdevColor);
% hold on

spikeInfo = {};

start_time = 1;
end_time = size(timeVector);
start_of_spike = 1;
end_of_spike = 1;
potential_end = 1;
inMiniSpike = 0;
inBigSpike = 0;
spikeFreqThresh = .005 * sr; %seconds threshhold * sample rate
timeSinceSpike = spikeFreqThresh + 1; % starting large enough to make the first spike count

%below finds spikes that are ABOVE the threshold
for i = start_time:end_time
    if(y(i) > spikeThresh && inMiniSpike == 0 && inBigSpike == 0 && timeSinceSpike > spikeFreqThresh)
        start_of_spike = timeVector(i);
        timeSinceSpike = 0;
        inMiniSpike = 1;
        inBigSpike = 1;
        potential_end = 1;
    elseif(y(i) > spikeThresh && inMiniSpike == 0 && inBigSpike == 1 && timeSinceSpike < spikeFreqThresh)
        timeSinceSpike = 0;
        inMiniSpike = 1;
        potential_end = 1;
    elseif(y(i) < spikeThresh && inMiniSpike == 0 && inBigSpike == 1 && timeSinceSpike > spikeFreqThresh)
        end_of_spike = potential_end;
        inMiniSpike = 0;
        inBigSpike = 0;
        timeOfSpike = findAboveSpike(start_of_spike, end_of_spike, y);
        if(timeOfSpike == 0)
            fprintf("No spike found between %5.5f and %5.5f\n", start_of_spike, end_of_spike);
            %above line should never print, there should always be a spike found
        else
            lastScent = findLastScent(timeOfSpike, linaloolConcTrace, ethanolConcTrace);
            if(strcmp(lastScent(1,1),""))
                fprintf("ERROR NO LAST SCENT FOUND");
            end
            spikeInfo = [spikeInfo; {start_of_spike, end_of_spike,timeOfSpike/sr,v(timeOfSpike), lastScent{1,1}, lastScent{1,2}, lastScent{1,3}/sr}]; 
            %^^add a new row of spike data, contains time when spike passes
            %cutoff line at the start and end, the time of the local peak,
            %and the value at that time

        end
        start_of_spike = 1;
        end_of_spike = 1;
    elseif(y(i) < spikeThresh && inMiniSpike == 1 && inBigSpike == 1)
        potential_end = timeVector(i);
        inMiniSpike = 0;
        timeSinceSpike = 0;
    end
    timeSinceSpike = timeSinceSpike + 1;
end

%below finds spikes that are BELOW the threshold
numberOfTopSpikes = size(spikeInfo);
start_time = 1;
end_time = size(timeVector);
start_of_spike = 1;
end_of_spike = 1;
potential_end = 1;
inMiniSpike = 0;
inBigSpike = 0;
spikeFreqThresh = .005 * sr; %seconds threshhold * sample rate
timeSinceSpike = spikeFreqThresh + 1; % starting large enough to make the first spike count
for i = start_time:end_time
    if(y(i) < (-1)*spikeThresh && inMiniSpike == 0 && inBigSpike == 0 && timeSinceSpike > spikeFreqThresh)
        start_of_spike = timeVector(i);
        timeSinceSpike = 0;
        inMiniSpike = 1;
        inBigSpike = 1;
        potential_end = 1;
    elseif(y(i) < (-1)*spikeThresh && inMiniSpike == 0 && inBigSpike == 1 && timeSinceSpike < spikeFreqThresh)
        timeSinceSpike = 0;
        inMiniSpike = 1;
        potential_end = 1;
    elseif(y(i) > (-1)*spikeThresh && inMiniSpike == 0 && inBigSpike == 1 && timeSinceSpike > spikeFreqThresh)
        end_of_spike = potential_end;
        inMiniSpike = 0;
        inBigSpike = 0;
        timeOfSpike = findBelowSpike(start_of_spike, end_of_spike, y);
        if(timeOfSpike == 0)
            fprintf("No spike found between %5.5f and %5.5f\n", start_of_spike, end_of_spike);
            %above line should never print, there should always be a spike found
        else
            lastScent = findLastScent(timeOfSpike, linaloolConcTrace, ethanolConcTrace);
            if(strcmp(lastScent(1,1),""))
                fprintf("ERROR NO LAST SCENT FOUND");
            end
            spikeInfo = [spikeInfo; {start_of_spike, end_of_spike,timeOfSpike/sr,v(timeOfSpike), lastScent{1,1}, lastScent{1,2}, lastScent{1,3}/sr}]; 
            %^^add a new row of spike data, contains time when spike passes
            %cutoff line at the start and end, the time of the local peak,
            %and the value at that time

        end
        start_of_spike = 1;
        end_of_spike = 1;
    elseif(y(i) > (-1)*spikeThresh && inMiniSpike == 1 && inBigSpike == 1)
        potential_end = timeVector(i);
        inMiniSpike = 0;
        timeSinceSpike = 0;
    end
    timeSinceSpike = timeSinceSpike + 1;
end

numberOfBottomSpikes = size(spikeInfo) - numberOfTopSpikes;

%marks the spikes on the plot of the voltages
%red for linalool, red for ethanool
numSpikes = size(spikeInfo);
numSpikes = numSpikes(1,1);
for i = 1:numSpikes
   % markSize = spikeInfo{i,7}/(31/6);
   % marSize = fix(markSize);
    if(strcmp(spikeInfo{i,5},"linalool"))
        plot(spikeInfo{i,3}, spikeInfo{i,4}, 'r*','MarkerSize',3);
    elseif(strcmp(spikeInfo{i,5},"ethanol"))
        plot(spikeInfo{i,3}, spikeInfo{i,4}, 'b*','MarkerSize',3);
    end
    hold on;
end

newFile = strcat(beeOutputFile, '\fullGraph');

print(newFile, '-dpng');

%% See periodogram - don't think this works yet
periodogram(vFilt, sr);

newFile = strcat(beeOutputFile, '\fft');

print(newFile, '-dpng');

%%
%this script finds the most recent scent from the spike, and makes a
%histogram of how long it has been since the scent presentation
linaloolScentPowers = {};
ethanolScentPowers = {};

for i = 1:numSpikes
   if(strcmp(spikeInfo{i,5},"linalool"))
        linaloolScentPowers = [linaloolScentPowers; {spikeInfo{i,7}}];
   elseif(strcmp(spikeInfo{i,5},"ethanol"))
        ethanolScentPowers = [ethanolScentPowers ; {spikeInfo{i,7}}];
   end
end

linaloolScentPowers = cell2mat(linaloolScentPowers);
ethanolScentPowers = cell2mat(ethanolScentPowers);

overHistcTrans(linaloolScentPowers,ethanolScentPowers,1);

newFile = strcat(beeOutputFile, '\histogram');

print(newFile, '-dpng');

%%
%this script will show the spikes overlayed atop eachother on a plot
figure;

for i = 1:size(spikeInfo) %Below, linalool, low
    if(spikeInfo{i,4} > 0 && spikeInfo{i,5} == 'linalool' && spikeInfo{i,6} == 11)
        spikeMidIndex = spikeInfo{i,3} * sr;
        spikeBounds = 1 * sr;

        spikeStartIndex = spikeMidIndex - spikeBounds;
        if(spikeStartIndex < 0)
            spikeStartIndex = 0;
        end
        spikeEndIndex = spikeMidIndex + spikeBounds;
        lengthOfSpike = spikeEndIndex - spikeStartIndex;
        
        additionalLength = size(v(spikeStartIndex:spikeEndIndex)) - size(timeVector(1:lengthOfSpike));
        lengthOfSpike = lengthOfSpike + additionalLength;
        
        overlay = plot(timeVector(1:lengthOfSpike), v(spikeStartIndex:spikeEndIndex), 'r');
        overlay.Color(4) = 0.05;
        hold on;
    end
end
xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Low Linalool Above Spikes");

newFile = strcat(beeOutputFile, '\lowTopLinaloolTrace');

print(newFile, '-dpng');

figure;

for i = 1:size(spikeInfo) %Below, linalool, med
    if(spikeInfo{i,4} > 0 && spikeInfo{i,5} == 'linalool' && spikeInfo{i,6} == 33)
        spikeMidIndex = spikeInfo{i,3} * sr;
        spikeBounds = 1 * sr;

        spikeStartIndex = spikeMidIndex - spikeBounds;
        spikeEndIndex = spikeMidIndex + spikeBounds;
        lengthOfSpike = spikeEndIndex - spikeStartIndex;
        
        additionalLength = size(v(spikeStartIndex:spikeEndIndex)) - size(timeVector(1:lengthOfSpike));
        lengthOfSpike = lengthOfSpike + additionalLength;
        
        overlay = plot(timeVector(1:lengthOfSpike), v(spikeStartIndex:spikeEndIndex), 'r');
        overlay.Color(4) = 0.05;
        hold on;
    end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Medium Linalool Above Spikes");

newFile = strcat(beeOutputFile, '\mediumTopLinaloolTrace');

print(newFile, '-dpng');

figure;

for i = 1:size(spikeInfo) %Below, linalool, high
    if(spikeInfo{i,4} > 0 && spikeInfo{i,5} == 'linalool' && spikeInfo{i,6} == 100)
        spikeMidIndex = spikeInfo{i,3} * sr;
        spikeBounds = 1 * sr;
        startDisplacement = 1;
       
        spikeStartIndex = spikeMidIndex - spikeBounds;
        if(spikeStartIndex < 1)
            startDisplacement = -1 * spikeStartIndex;
            spikeStartIndex = 1;
        end
        spikeEndIndex = spikeMidIndex + spikeBounds;
        lengthOfSpike = spikeEndIndex - spikeStartIndex;
        
        additionalLength = size(v(spikeStartIndex:spikeEndIndex)) - size(timeVector(1:lengthOfSpike));
        lengthOfSpike = lengthOfSpike + additionalLength;
        
        overlay = plot(timeVector(startDisplacement:lengthOfSpike), v(spikeStartIndex:spikeEndIndex), 'r');
        overlay.Color(4) = 0.05;
        hold on;
    end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("High Linalool Above Spikes");

newFile = strcat(beeOutputFile, '\highTopLinaloolTrace');

print(newFile, '-dpng');

%this script will show the spikes overlayed atop eachother on a plot
figure;

for i = 1:size(spikeInfo) %Below, ethanol, unfiltered
    if(spikeInfo{i,4} > 0 && spikeInfo{i,5} == 'ethanol')
        spikeMidIndex = spikeInfo{i,3} * sr;
        spikeBounds = 1 * sr;

        spikeStartIndex = spikeMidIndex - spikeBounds;
        spikeEndIndex = spikeMidIndex + spikeBounds;
        lengthOfSpike = spikeEndIndex - spikeStartIndex;
        
        additionalLength = size(v(spikeStartIndex:spikeEndIndex)) - size(timeVector(1:lengthOfSpike));
        lengthOfSpike = lengthOfSpike + additionalLength;
        
        overlay = plot(timeVector(1:lengthOfSpike), v(spikeStartIndex:spikeEndIndex), 'b');
        overlay.Color(4) = 0.05;
        hold on;
    end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Ethanol Above Spikes");

newFile = strcat(beeOutputFile, '\allTopEthanolTrace');

print(newFile, '-dpng');

%%

figure;

for i = 1:size(linaloolMat)%low linalool window
   if(linaloolMat(i,3) ==  11)
       scentStart = linaloolMat(i,1) + 1; %adding 1 because vectors start at 1
       scentEnd = linaloolMat(i,2) + 1;
       overlay = plot(timeVector(1:1001), v(scentStart:scentEnd), 'r');
       overlay.Color(4) = 0.2;
       hold on;
   end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Low Linalool Scent Window");

newFile = strcat(beeOutputFile, '\lowLinaloolWindow');

print(newFile, '-dpng');

figure;

for i = 1:size(linaloolMat)%med linalool window
   if(linaloolMat(i,3) ==  33)
       scentStart = linaloolMat(i,1) + 1; %adding 1 because vectors start at 1
       scentEnd = linaloolMat(i,2) + 1;
       overlay = plot(timeVector(1:1001), v(scentStart:scentEnd), 'r');
       overlay.Color(4) = 0.2;
       hold on;
   end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Medium Linalool Scent Window");

newFile = strcat(beeOutputFile, '\medLinaloolWindow');

print(newFile, '-dpng');


figure;

for i = 1:size(linaloolMat)%high linalool window
   if(linaloolMat(i,3) ==  100)
       scentStart = linaloolMat(i,1) + 1; %adding 1 because vectors start at 1
       scentEnd = linaloolMat(i,2) + 1;
       overlay = plot(timeVector(1:1001), v(scentStart:scentEnd), 'r');
       overlay.Color(4) = 0.2;
       hold on;
   end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("High Linalool Scent Window");

newFile = strcat(beeOutputFile, '\highLinaloolWindow');

print(newFile, '-dpng');


figure;

for i = 1:size(ethanolMat)%ethanol window
   if(ethanolMat(i,3) ==  100)
       scentStart = ethanolMat(i,1) + 1; %adding 1 because vectors start at 1
       scentEnd = ethanolMat(i,2) + 1;
       overlay = plot(timeVector(1:1001), v(scentStart:scentEnd), 'b');
       overlay.Color(4) = 0.2;
       hold on;
   end
end

xlabel("Relative Time (s)");
ylabel("Responses (uV)");
title("Ethanol Scent Window");

newFile = strcat(beeOutputFile, '\ethanolWindow');

print(newFile, '-dpng');

close all
end
function timeOfSpike = findAboveSpike(start_point, end_point, vFilt)
start_point = fix(start_point * 20000); %makes it a integer so it can be used in a loop
end_point = fix(end_point * 20000);
%finds the top of the spike
for i = start_point:end_point
    peak_top = true;
    for j = i:end_point
        if(vFilt(i) < vFilt(j))
            peak_top = false;
            break;
        end
    end
    if(peak_top == true)
        timeOfSpike = i;
        return;
    end
end
timeOfSpike = 0;
return
end

function timeOfSpike = findBelowSpike(start_point, end_point, vFilt)
start_point = fix(start_point * 20000); %makes it a integer so it can be used in a loop
end_point = fix(end_point * 20000);
%fprintf(" Start: %lu \n End: %lu \n", start_point, end_point);
%finds the top of the spike
for i = start_point:end_point
    peak_top = true;
    for j = i:end_point
        if(vFilt(i) > vFilt(j))
            peak_top = false;
            break;
        end
    end
    if(peak_top == true)
        timeOfSpike = i;
        return;
    end
end
timeOfSpike = 0;
return
end

function lastScent = findLastScent(timeOfSpike, linaloolConcTrace, ethanolConcTrace)
for j = timeOfSpike:-1:1
    if(linaloolConcTrace(j) ~= 0)
        lastScent = {"linalool", linaloolConcTrace(j), timeOfSpike - j};
        return;
    elseif(ethanolConcTrace(j) ~= 0)
        lastScent = {"ethanol", ethanolConcTrace(j), timeOfSpike - j};
        return;
    end
end
return;
end

function [periodogramH] = periodogram(inputVector, sr)

% Plot periodogram

% N = length(inputVector);
% xdft = fft(inputVector);
% xdft = xdft(1:N/2+1);
% psdx = (1/(sr*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:sr/length(inputVector):sr/2;
% 
% periodogramH = figure;
% plot(freq,10*log10(psdx))
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

N = length(inputVector);
Y = fft(inputVector);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = sr*(0:(N/2))/N;
%display(f);
periodogramH = figure;
plot(f,P1)
hold on
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
xlim([0 100])
end

function [odorConcTrace] = makeOdorTrace(tVector, sr, odorMat, varargin)
% Creates a vector of zeros of the same size as timeVector, which must be
% passed in as tVector.  Changes the zeros to match the odor concentrations
% during odor delivery. Outputs a trace. Can optionally add shaded areas to
% the summary plot.

% Example: linaloolConcTrace = makeOdorTrace(timeVector, sr, linaloolMat, 'shadeColor', [255 0 0]);
%

% Optional arguments
shadeColor = [];

pvpmod(varargin)

odorConcTrace = zeros(size(tVector));
% Convert from ms to samples
odorMat(:,1:2) = odorMat(:,1:2).*sr/1000 + 1; %must add 1 because time 0 starts at sample 1
for i = 1:size(odorMat, 1)
    odorConcTrace(odorMat(i,1):odorMat(i,2)) = odorMat(i,3);
end
if isempty(shadeColor)
    %do nothing
else
    for i = 1:size(odorMat,1)
        AlphaVal = .6;
        
        x = [tVector(odorMat(i,1)) tVector(odorMat(i,2)) ];
        posy = [500 500].*(odorMat(i,3)/100);
        negy = [-500 -500].*(odorMat(i,3)/100);
        shadeH(i) = area(x,posy);
        shadeH(i).FaceColor = shadeColor;
        shadeH(i).FaceAlpha = AlphaVal;
        shadeH(i).EdgeAlpha = 0;
        shadeH(i) = area(x,negy);
        shadeH(i).FaceColor = shadeColor;
        shadeH(i).FaceAlpha = AlphaVal;
        shadeH(i).EdgeAlpha = 0;
    end
end

end
