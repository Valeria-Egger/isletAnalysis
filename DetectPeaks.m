function [signals, averageFilter, stdFilter] = DetectPeaks(y, lag, threshold, influence, Minimum_signal)

signals = zeros(length(y), 1);
filteredY = zeros(size(y));
filteredY(1:lag+1) = y(1:lag+1);
averageFilter(lag+1, 1) = mean(y(1:lag+1));
stdFilter(lag+1, 1) = std(y(1:lag+1));

for i=lag+2:length(y)
    if y(i)<Minimum_signal
        signals(i) = 0;
    else 
        if abs(y(i)-averageFilter(i-1))>threshold*stdFilter(i-1)
            if y(i)>averageFilter(i-1)
                signals(i) = 1;
            else
            signals(i) = 0;
            end
            filteredY(i) = influence*y(i)+(1-influence)*filteredY(i-1);
        else
            signals(i) = 0; 
            filteredY(i) = y(i);
        end 
    end
    averageFilter(i, 1) = mean(filteredY(i-lag:i));
    stdFilter(i, 1) = std(filteredY(i-lag:i));
end
end
