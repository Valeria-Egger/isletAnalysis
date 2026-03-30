% Load min-max scaled Data

[file, path] = uigetfile("*.csv", "Select a csv file");
[~, name, ~] = fileparts(file);
dataset = readtable(fullfile(path, file));

%Preallocation
numCols = width(dataset);
AllPeakStarts = cell(numCols, 1);
AllPeakProminence = cell(numCols, 1);
AllPeakHeights = cell(numCols, 1);
AllPeakWidths = cell(numCols, 1);
AllPeakEven = cell(numCols, 1);
AllPeakSharpSure = cell(numCols, 1);
AllInflectionsRise = cell(numCols, 1);
AllInflectionsDecay = cell(numCols, 1);
SymmetryMatrix = cell(numCols, 1);
AllAUCs = cell(numCols, 1);
allSignals = cell(1, numCols);
AllCounts = cell(1, numCols);
PeakCounts = zeros(1, numCols);
AllInflectionPointRise = cell(1, numCols);
AllInflectionPointDecay = cell(1, numCols);
AllPeakIndex = cell(1, numCols);
AllRawSignal = cell(1, numCols);
CoherenceMatrix = zeros(numCols, numCols);
LeadershipScore = zeros(numCols, numCols);

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
for i = 2:numCols
    rawSignal = dataset{:, i};
    rawSignal = rawSignal(:);
    AllRawSignal{i} = rawSignal;
    [signals, avg, dev] = DetectPeaks(rawSignal, lag, threshold, influence, Minimum_signal);
    peaks = extractPeaks(rawSignal, signals);

    AllPeakStarts{i} = peaks.Start; 
    AllPeakProminence{i} = peaks.Prominence;
    AllPeakWidths{i} = peaks.Widths;
    AllPeakHeights{i} = peaks.Heights; 
    AllPeakEven{i} = peaks.even;
    AllPeakSharpSure{i} = peaks.Sharp;
    AllInflectionsRise{i} = peaks.riseInflectionIndex;
    AllInflectionsDecay{i} = peaks.decayInflectionIndex;
    SymmetryMatrix{i} = AllInflectionsRise{i} - AllInflectionsDecay{i};
    AUC = AUCPerWindow(rawSignal, 180);
    AllAUCs{i} = AUC';

    AllInflectionPointRise{i} = peaks.riseInflectionPoint;
    AllInflectionPointDecay{i} = peaks.decayInflectionPoint;
    AllPeakIndex{i} = peaks.PeakIndex;
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

Wavelet{exampleCell} = ComputeWavelet(AllRawSignal{exampleCell});
PlotWavelet(Wavelet{exampleCell}, exampleCell);

P = cell(1, numCols);
W = cell(1, numCols);



%N = numel(W);
for i = 2:numCols
    W{i} = ComputeWavelet(AllRawSignal{i});
    P{i} = abs(W{i}.cfs).^2;
end
Pstack = cat(3, P{:});
Pmean = mean(Pstack, 3);
Pfreq = mean(Pmean, 2);
numTimes = size(Pmean, 2);
time = (0:numTimes-1) * 1;
freqs = W{2}.freqs;

for i = 2:numCols
    for j = 2:numCols
    coh_threshold = 0.4;
    x = AllRawSignal{i};
    z = AllRawSignal{j};
    [wcoh, wcs, period, coi] = wcoherence(x, z);
    period = period(:);
    coi = coi(:)';
    mask = wcoh;
    mask(coi < period) = NaN;
    meanCoh = mean(wcoh(:), 'omitnan');
    CoherenceMatrix(i, j) = meanCoh;

    phase = angle(wcs);
    valid = (wcoh > coh_threshold) & (coi >= period);
    validPhase = phase(valid);
    xLeads = sum(validPhase > 0);
    yLeads = sum(validPhase < 0);
    total = xLeads + yLeads;
    if total == 0
        score = NaN;
    else
        score = (xLeads-yLeads)/total;
    end
    LeadershipScore(i, j) = score;
    end
end

%Plotting
figure;
imagesc(CoherenceMatrix(2:end, 2:end));
colormap("turbo");
colorbar;
xlabel("i cell");
ylabel("j cell");
title("Mean wavelet coherence strength between cell i and j");

figure;
imagesc(LeadershipScore(2:end, 2:end));
colormap("turbo");
colorbar;
xlabel("i cell");
ylabel("j cell");
title("Leadershipscore between cell i and j");

figure;
imagesc(time, freqs, Pmean);
axis xy;
colormap("turbo");
colorbar;
xlabel("time");
ylabel("Frequency");
title("population Mean wavelet Power");

figure;
plot(freqs, Pfreq, 'LineWidth', 2);
xlabel("Frequency");
ylabel("Mean Power");
title("Population Power spectrum");
grid on;

[acf, lags] = ComputeAutocorr(allSignals{exampleCell});
pos = lags >= 0;
figure;
    plot(lags(pos), acf(pos), 'LineWidth', 1.5);
    xlabel('Lag');
    ylabel('Autocorrelation');
    title(['Autocorrelation (positive lags) for Cell', num2str(exampleCell)])


raw = dataset{:, exampleCell};
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

AllPeakHeights = cellfun(@(x) x(:), AllPeakHeights, 'UniformOutput', false);
allHeights = vertcat(AllPeakHeights{:});
edges = 0:0.2:max(allHeights);
countsH = histcounts(allHeights, edges);

figure;
bar(edges(1:end-1), countsH, 'histc');
xlabel("Peak Heights");
ylabel("Count");
title("Distribution of Peak Heights")

AllPeakWidths = cellfun(@(x) x(:), AllPeakWidths, 'UniformOutput', false);
allWidths = vertcat(AllPeakWidths{:});
edges = 0:10:max(allWidths);
countsW = histcounts(allWidths, edges);

figure;
bar(edges(1:end-1), countsW, 'histc');
xlabel("Peak Widths");
ylabel("Count");
title("Distribution of Peak Widths")

nonEmpty = ~cellfun(@isempty, AllCounts);
Peak_Matrix_count = cell2mat(AllCounts(nonEmpty)');
total_per_window = sum(Peak_Matrix_count, 1);
figure;
bar(total_per_window);
title('total peaks per time window');

nonEmpty = ~cellfun(@isempty, AllAUCs);
AUC_Matrix_count = cell2mat(AllAUCs(nonEmpty));
total_AUC_per_window = sum(AUC_Matrix_count, 1);
figure;
bar(total_AUC_per_window);
title('total AUC per time window');

figure;
bar(PeakCounts);
title('peaks per cell');

AllPeakProminence = cellfun(@(x) x(:), AllPeakProminence, 'UniformOutput', false);
allProminence = vertcat(AllPeakProminence{:});
edges = 0:0.2:max(allProminence);
countsP = histcounts(allProminence, edges);
figure;
bar(edges(1:end-1), countsP, 'histc');
xlabel("Peak Prominence");
ylabel("Count");
title("Distribution of Peak Prominence")

counts = zeros(numCols, 2);
for i = 1:numCols
    vector = AllPeakEven{i};
    counts(i, 1) = sum(vector);
    counts(i, 2) = sum(~vector);
end

figure;
bar(counts);
set(gca, 'XTick', 2:numCols);
legend({'Symmetric', 'Asymmetric'});
ylabel('Number of Peaks');
title("Symmetry over the dataset")

countsSharp = zeros(numCols, 2);
vector2 = zeros(numCols, 2);



for i = 1:numCols
    vector2 = AllPeakSharpSure{i};
    countsSharp(i, 1) = sum(vector2);
    countsSharp(i, 2) = sum(~vector2);
end

figure;
bar(countsSharp, "stacked");
set(gca, 'XTick', 1:numCols);
legend({'Sharp', 'The opposite of Sharp'});
ylabel('Number of Peaks');
title("Sharpness over the dataset")



maxLen = max(cellfun(@length, SymmetryMatrix));
SymPad = nan(numCols, maxLen);

for i = 1:numCols
    L = length(SymmetryMatrix{i});
    SymPad(i,1:L) = SymmetryMatrix{i};
end

figure;
boxplot(SymPad', 'Labels', 1:numCols);
xlabel('Cell #');
ylabel('Symmetry Score');
title('Symmetry Distribution per Cell');

%symmetry score >0 longer rise than decay, peak leans right
%symmetry score <0 longer decay, peak leans left
%symmetry score = 0 rise and decay equal
%we are back with unnecessary documentation somewhere

RiseInflectionPointEx = AllInflectionPointRise{exampleCell};
DecayInflectionPointEx = AllInflectionPointDecay{exampleCell};


figure;
hold on;
plot(raw, 'k', 'LineWidth', 1.2);
scatter(AllPeakIndex{exampleCell}, raw(AllPeakIndex{exampleCell}), 60, 'r', 'filled', ...
    'DisplayName', 'PeakMaxima');
validRise = ~isnan(RiseInflectionPointEx);
validDecay = ~isnan(DecayInflectionPointEx);
scatter(RiseInflectionPointEx(validRise), raw(RiseInflectionPointEx(validRise)),...
    50, 'b', 'filled', 'DisplayName', 'RiseInflection');
scatter(DecayInflectionPointEx(validDecay), raw(DecayInflectionPointEx(validDecay)),...
    50, 'green', 'filled', 'DisplayName', 'DecayInflection');
legend show;
title('Peak Detection and Inflection Points');
xlabel('Sample Index');
ylabel('Signal Amplitude');
hold off;


%not a great fan of that plot to be honest but however
RiseInflectionIndexEx = AllInflectionsRise{exampleCell};
DecayInflectionIndexEx = AllInflectionsDecay{exampleCell};
figure; hold on;
scatter(RiseInflectionIndexEx, DecayInflectionIndexEx, 60, 'filled');
plot([0 max(RiseInflectionIndexEx)], [0 max(DecayInflectionIndexEx)], 'k--'); % symmetry line
xlabel('Rise-side distance');
ylabel('Decay-side distance');
title('Symmetry Scatter Plot');
axis equal;
grid on;
hold off;

%save the data
excelFile = name + "AnalysisResults.xlsx";
CellIDs = "Cell" + (1:numCols)';
PeakCounts_table = table(CellIDs, PeakCounts(:), ...
    'VariableNames', {'CellsIDs', 'PeakCount'});
writetable(PeakCounts_table, excelFile, 'Sheet', 'PeakCount');

rows = [];


for i = 1:numCols
    heights = AllPeakHeights{i};
    widths = AllPeakWidths{i};
    prominence = AllPeakProminence{i};

    numPeaks = length(heights);
    CellID = repmat(i, numPeaks, 1);
    PeakIndex = (1:numPeaks)';
    rows = [rows; table(CellID, PeakIndex, heights(:), widths(:), prominence(:),...
        'VariableNames', {'CellID', 'CellIndex', 'CellHeight', 'CellWidth', 'CellProminence'})];
end
writetable(rows, excelFile, 'Sheet', 'PeakDetails');

%numCells = size(AUC_Matrix_count, 1);
%CellIDs = "Cell" + (1:numCells)';

%AUC_table = array2table(AUC_Matrix_count, ...
%    'RowNames', CellIDs);

%writetable(AUC_table, excelFile, 'Sheet', 'AUC_per_window', 'WriteRowNames', true);


%numCells = size(Peak_Matrix_count, 1);
%CellIDs = "Cell" + (1:numCells)';

%Peak_table = array2table(Peak_Matrix_count, ...
%    'RowNames', CellIDs);

%writetable(Peak_table, excelFile, 'Sheet', 'Peak_per_window', 'WriteRowNames', true);


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







