% Load min-max scaled Data

[file, path] = uigetfile("*.csv", "Select a csv file");
dataset = readtable(fullfile(path, file));

AllPeakStarts = cell(numCols, 1);
AllPeakProminence = cell(numCols, 1);
AllPeakHeights = cell(numCols, 1);
AllPeakWidths = cell(numCols, 1);

numCols = width(dataset);
allSignals = cell(1, numCols);
AllCounts = cell(1, numCols);
PeakCounts = zeros(1, numCols);

% Settings
lag = 100;
threshold = 2;
influence = 0.01;
Minimum_signal = 0.2;
exampleCell = 7;

%Calculate the number of peaks detected
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

%Peak feature extraction
for i = 1:numCols
    rawSignal = dataset{:, i};
    [signals, avg, dev] = DetectPeaks(rawSignal, lag, threshold, influence, Minimum_signal);
    [PeakStarts, PeakEnds, PeakHeights, PeakWidths, peakProminence] = extractPeaks(rawSignal, signals);

    AllPeakStarts{i} = PeakStarts; 
    AllPeakProminence{i} = peakProminence(:);
    AllPeakWidths{i} = PeakWidths(:);
    AllPeakHeights{i} = PeakHeights(:); 
end

%define time windows (this may have to be adjusted for your time windows)
%for i = 2:numCols
%winSize = 180;
%signals = allSignals{i}(:);
%AllCounts{i} = countPeaksPerWindow(signals, winSize);
%end

%disp(AllCounts);

%Calculate the number of peaks per islet (so for the whole dataset)
%right now just one dataset because I am tired

totalPeaks = sum(PeakCounts);
fprintf('Total number of peaks detected across all columns: %d\n', totalPeaks);

for i = 2:numCols
fprintf('Number of peaks detected in column %d: %d\n', i, PeakCounts(i));
end

%Wavelet and autocorrelation for oscillation detection
W{exampleCell} = ComputeWavelet(allSignals{exampleCell});
PlotWavelet(W{exampleCell}, exampleCell);

%Plotting
[acf, lags] = ComputeAutocorr(allSignals{exampleCell});
pos = lags >= 0;
figure;
    plot(lags(pos), acf(pos), 'LineWidth', 1.5);
    xlabel('Lag');
    ylabel('Autocorrelation');
    title(['Autocorrelation (positive lags) for Cell', num2str(exampleCell)])


raw = dataset{:, 7};
[signals1,avg1,dev1] = DetectPeaks(raw,lag,threshold,influence, Minimum_signal);

figure; subplot(2,1,1); hold on;
x = 1:length(raw); 
t = lag+1:length(raw);
area(x(t),avg1(t)+threshold*dev1(t),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
area(x(t),avg1(t)-threshold*dev1(t),'FaceColor',[1 1 1],'EdgeColor','none');
plot(x(t),avg1(t),'LineWidth',1,'Color','cyan','LineWidth',1.5);
plot(x(t),avg1(t)+threshold*dev1(t),'LineWidth',1,'Color','green','LineWidth',1.5);
plot(x(t),avg1(t)-threshold*dev1(t),'LineWidth',1,'Color','green','LineWidth',1.5);
plot(1:length(raw),raw,'b');
subplot(2,1,2);
stairs(signals1,'r','LineWidth',1.5); ylim([-1.5 1.5]);

allHeights = vertcat(AllPeakHeights{:});
edges = 0:0.2:max(allHeights);
countsH = histcounts(allHeights, edges);

figure;
bar(edges(1:end-1), countsH, 'histc');
xlabel("Peak Heights");
ylabel("Count");
title("Distribution of Peak Heights")

allWidths = vertcat(AllPeakWidths{:});
edges = 0:10:max(allWidths);
countsW = histcounts(allWidths, edges);

figure;
bar(edges(1:end-1), countsW, 'histc');
xlabel("Peak Widths");
ylabel("Count");
title("Distribution of Peak Widths")

%nonEmpty = ~cellfun(@isempty, AllCounts);
%Matrix_count = cell2mat(AllCounts(nonEmpty)');
%total_per_window = sum(Matrix_count, 1);
%figure;
%bar(total_per_window);
%title('total peaks per time window');

figure;
histogram(PeakCounts);
title('peaks per cell');


allProminence = vertcat(AllPeakProminence{:});
edges = 0:0.2:max(allProminence);
countsP = histcounts(allProminence, edges);
figure;
bar(edges(1:end-1), countsP, 'histc');
xlabel("Peak Prominence");
ylabel("Count");
title("Distribution of Peak Prominence")


%figure
%for i = 1:numCols 
%    if isempty(AllPeakStarts{i}) 
%        fprintf("Skipping Cell %d (no peaks)\n", i); 
%        continue 
%    end 
%    subplot(5, 4, i); 
%    stem(AllPeakStarts{i}, AllPeakProminence{i}, 'filled'); 
%    xlabel('Sample index'); 
%    ylabel('Prominence'); 
%    title(['Peak Prominence for Cell ' num2str(i)]); 
%end




