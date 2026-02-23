function [acf, lags] = ComputeAutocorr(y)
    y = y-mean(y);
    [acf, lags] = xcorr(y, 'coeff');
    pos = lags>=0;

    figure;
    plot(lags(pos), acf(pos), 'LineWidth', 1.5);
    xlabel('Lag');
    ylabel('Autocorrelation');
    title(['Autocorrelation (positive lags) for Cell', num2str(i)])
end