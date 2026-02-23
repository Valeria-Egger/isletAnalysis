% Source - https://stackoverflow.com/a/54507329
% Posted by Jean-Paul, modified by community. See post 'Timeline' for change history
% Retrieved 2026-02-12, License - CC BY-SA 4.0

% Data

dataset = readtable("scaled/Scaled_Results_ablated_islet3.csv");

numCols = width(dataset);
allSignals = cell(1, numCols);
AllCounts = cell(1, numCols);
PeakCounts = zeros(1, numCols);
% Settings
lag = 100;
threshold = 2;
influence = 0.01;
Minimum_signal = 0.2;

 %Calculate the number of peaks detected for the current signal
 %I will later rewrite this one with a rising edge detector as I will do
 %below
for i = 1:numCols
    y = dataset{:, i};
    [signals,avg,dev] = DetectPeaks(y,lag,threshold,influence, Minimum_signal);
    allSignals{i} = signals;
    numPeaks = 0;
    for k = 2:length(signals)
    if signals(k)==1 && signals(k-1) == 0
        numPeaks = numPeaks +1;
    end 
    end
    PeakCounts(i) = numPeaks;
end

for i = 1:numCols
    rawSignal = dataset{:, i};
    [signals, avg, dev] = DetectPeaks(rawSignal, lag, threshold, influence, Minimum_signal);
    allSignals{i} = signals;

    k = 1;
    PeakStarts = [];
    PeakEnds = [];
    PeakHeights = [];
    while k < length(signals)
        if signals(k) == 1
            startIdx = k;
            while k <= length(signals) && signals(k) == 1
                k = k +1;
            end
            endIdx = k-1;
            PeakStarts(end+1) = startIdx;
            PeakEnds(end+1) = endIdx;

            PeakHeights(end+1) = max(rawSignal(startIdx:endIdx));
            PeakWidths = PeakEnds - PeakStarts+1;

        end
        k = k+1;
    end

    peakProminence = zeros(1, length(PeakHeights));
    for p = 1:length(PeakHeights)
        s = PeakStarts(p);
        e = PeakEnds(p);

        left = min(rawSignal(1:s));
        right = min(rawSignal(e:end));

        peakProminence(p) = PeakHeights(p) - max(left, right);
    end
    PeakStarts_all{i} = PeakStarts; 
    PeakProminence_all{i} = peakProminence;
end

%for i = 1:numCols 
%    if isempty(PeakStarts_all{i}) 
%        fprintf("Skipping Cell %d (no peaks)\n", i); 
%        continue 
%    end 
%    figure; 
%    stem(PeakStarts_all{i}, PeakProminence_all{i}, 'filled'); 
%    xlabel('Sample index'); 
%    ylabel('Prominence'); 
%    title(['Peak Prominence for Cell ' num2str(i)]); 
%end

time = length(signals);
%adjust these as you needed if you have different conditions
for i = 2:numCols
winSize = 180;
signals = allSignals{i}(:);

N = floor(length(signals)/winSize)*winSize;
sigTrim = signals(1:N);

SignalWindow = reshape(sigTrim, winSize, []);

edges = SignalWindow(2:end, :) == 1 & SignalWindow(1:end-1, :) == 0; 
counts = sum(edges, 1);
AllCounts{i} = counts;
%disp(counts);
end
disp(AllCounts);

%Calculate the number of peaks per islet (so for the whole dataset)
%right now just one dataset because I am tired
totalPeaks = sum(PeakCounts);
fprintf('Total number of peaks detected across all columns: %d\n', totalPeaks);

for i = 2:numCols
fprintf('Number of peaks detected in column %d: %d\n', i, PeakCounts(i));
end

%total peak count per time window
nonEmpty = ~cellfun(@isempty, AllCounts);
Matrix_count = cell2mat(AllCounts(nonEmpty)');
total_per_window = sum(Matrix_count, 1);

%Wavelets = cell(1, numCols);
%for i = 2:numCols
%    y = dataset{:, i};
%    [cfs, freqs] = cwt(y, 'amor');
%    Wavelets{i}.cfs = cfs;
%    Wavelets{i}.freqs = freqs;
%end

%for i = 2:numCols
%   y = dataset{:, i};
%    figure;
%    cwt(y, 'amor');
%    title(['Wavelet Transform for Cell ' num2str(i)]);
%end

%try wavelet and autocorrelation here
W{i} = ComputeWavelet(y);
PlotWavelet(W{i}, 5);
ComputeAutocorr(y);




%iterate over all the cells then calculate number of peaks per cell, number
%of total peaks and number of peaks per time window per cell and per total
%peak detection
%number of peaks
%height of peaks
%length of peaks
%comparison wavelet and autocorrelation
%adjust minimal signal
%oscillation logic - autocorrelation may be better?
raw = dataset{:, 5};
[signals1,avg1,dev1] = DetectPeaks(raw,lag,threshold,influence, Minimum_signal);

%figure; subplot(2,1,1); hold on;
%x = 1:length(raw); 
%t = lag+1:length(raw);
%area(x(t),avg1(t)+threshold*dev1(t),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
%area(x(t),avg1(t)-threshold*dev1(t),'FaceColor',[1 1 1],'EdgeColor','none');
%plot(x(t),avg1(t),'LineWidth',1,'Color','cyan','LineWidth',1.5);
%plot(x(t),avg1(t)+threshold*dev1(t),'LineWidth',1,'Color','green','LineWidth',1.5);
%plot(x(t),avg1(t)-threshold*dev1(t),'LineWidth',1,'Color','green','LineWidth',1.5);
%plot(1:length(raw),raw,'b');
%subplot(2,1,2);
%stairs(signals1,'r','LineWidth',1.5); ylim([-1.5 1.5]);

%Histogram Peak Counts/per window
%figure('Color','w');   % clean new figure
%clf;                   % clear it just in case

%histogram(allCountsflat, 'BinMethod', 'integers', ...
%    'FaceColor',[0.2 0.4 0.8], 'EdgeColor','black');

%xlim([min(allCountsflat)-1, max(allCountsflat)+1]);
% title('Peak Counts per Window');
%xlabel('Peak Count');
%ylabel('Frequency');

%barplot to see total peaks per time window
figure;
bar(total_per_window);
title('total peaks per time window');

%histogram to see peaks per cell
figure;
histogram(PeakCounts);
title('peaks per cell');



