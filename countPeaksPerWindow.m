function counts = countPeaksPerWindow(signals, winSize)
N = floor(length(signals)/winSize)*winSize;
sigTrim = signals(1:N);

SignalWindow = reshape(sigTrim, winSize, []);

edges = SignalWindow(2:end, :) == 1 & SignalWindow(1:end-1, :) == 0; 
counts = sum(edges, 1);
end