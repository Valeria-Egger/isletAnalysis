function [acf, lags] = ComputeAutocorr(y)
    y = y-mean(y);
    [acf, lags] = xcorr(y, 'coeff');
end
