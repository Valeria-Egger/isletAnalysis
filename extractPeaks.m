function peaks = extractPeaks(rawSignal, signals)
   
k = 1;
    PeakStarts = [];
    PeakEnds = [];
    PeakHeights = [];
    PeakWidths = [];
    PeakIndex = [];

    while k < length(signals)
        if signals(k) == 1
            startIdx = k;
            while k <= length(signals) && signals(k) == 1
                k = k +1;
            end
            endIdx = k-1;
            PeakStarts(end+1) = startIdx;
            PeakEnds(end+1) = endIdx;
            [val, index] = max(rawSignal(startIdx:endIdx));
            PeakHeights(end+1) = val;
            PeakIndex(end + 1) = startIdx + index -1;
            PeakWidths(end+1) = endIdx - startIdx + 1;

        end
        k = k+1;
    end

    peakProminence = zeros(1, length(PeakHeights));
    AUC_left = zeros(1, length(PeakHeights));
    AUC_right = zeros(1, length(PeakHeights));
    PeakEven = zeros(1, length(PeakHeights));
    Symmetric = zeros(1, length(PeakHeights));
    Asymmetric = zeros(1, length(PeakHeights));
    
    for p = 1:length(PeakHeights)
        s = PeakStarts(p);
        e = PeakEnds(p);

        left = min(rawSignal(1:s));
        right = min(rawSignal(e:end));

        AUC_left(p) = trapz(rawSignal(s:PeakIndex(p)));
        AUC_right(p) = trapz(rawSignal(PeakIndex(p):e));

        %fprintf("AUC left is %d\n", AUC_left(p));
        %fprintf("AUC right is %d\n", AUC_right(p));
        DeltaAUC = abs(AUC_left(p) - AUC_right(p));
        PeakEven(p) = DeltaAUC < 2;

        Symmetric(p) = sum(PeakEven);
        Asymmetric(p) = ~sum(PeakEven);

        peakProminence(p) = PeakHeights(p) - max(left, right);
    end

%second derivative
secondDerivativ = zeros(size(rawSignal));
%approximation for second derivative for discrete numbers
secondDerivativ(2:end-1) = rawSignal(3:end) - 2*rawSignal(2:end-1) + rawSignal(1:end-2);

curvature = nan(1, length(PeakHeights));
rise_inflectionPoint = nan(1, length(PeakHeights));
decay_inflectionPoint = nan(1, length(PeakHeights));
rise_inflection_index = nan(1, length(PeakHeights));
decay_inflection_index = nan(1, length(PeakHeights));
PeakSharp = [];
   for p = 1:length(PeakHeights)
    s = PeakStarts(p);
    e = PeakEnds(p);
    PI = PeakIndex(p);
    curvature(p) = secondDerivativ(PI);
    %find inflection point at rise side (so from start to the actual
    %maximum
    rise_Index = s:PI;
    rise_SecDer = secondDerivativ(rise_Index);
    decay_Index = PI:e;
    decay_SecDer = secondDerivativ(decay_Index);
    signChange = find(diff(sign(rise_SecDer)) ~= 0, 1, 'first');
    if ~isempty(signChange)
        rise_inflectionPoint(p) = rise_Index(signChange);
        rise_inflection_index(p) = PI - rise_inflectionPoint(p);
    else
        rise_inflectionPoint(p) = NaN; 
        rise_inflection_index(p) = NaN;
    end
    signChange = find(diff(sign(decay_SecDer)) ~= 0, 1, 'first');
    if ~isempty(signChange)
        decay_inflectionPoint(p) = decay_Index(signChange);
        decay_inflection_index(p) = decay_inflectionPoint(p) - PI;
    else
        decay_inflectionPoint(p) = NaN;
        decay_inflection_index(p) = NaN;
    end
   end

SharpThreshold = prctile(curvature, 25);
PeakSharp = curvature < SharpThreshold;
PeakSharp = logical(PeakSharp)

peaks.Start = PeakStarts;
peaks.Ends = PeakEnds;
peaks.Heights = PeakHeights;
peaks.Widths = PeakWidths;
peaks.Prominence = peakProminence;
peaks.Sharp = PeakSharp;
peaks.decayInflectionPoint = decay_inflectionPoint;
peaks.decayInflectionIndex = decay_inflection_index;
peaks.riseInflectionPoint = rise_inflectionPoint;
peaks.riseInflectionIndex = rise_inflection_index;
peaks.curvature = curvature;
peaks.even = PeakEven;
peaks.AUCleft = AUC_left;
peaks.AUCright = AUC_right;
peaks.PeakIndex = PeakIndex;

end
