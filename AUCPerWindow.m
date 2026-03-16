function AUC = AUCPerWindow(y, winSize)
    N = floor(length(y)/winSize)*winSize;
    sigTrim = y(1:N);

    yWindow = reshape(sigTrim, winSize, []);
    numWindows = size(yWindow, 2)
    AUC = zeros(numWindows, 1);

    for i = 1:numWindows
        AUC(i) = trapz(yWindow(:, i));
    end
end
