function [PeakStarts, PeakEnds, PeakHeights, PeakWidths, peakProminence] = extractPeaks(rawSignal, signals)
   
k = 1;
    PeakStarts = [];
    PeakEnds = [];
    PeakHeights = [];
    PeakWidths = [];
    peakProminence = [];
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
            PeakWidths(end+1) = endIdx - startIdx + 1;

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
end