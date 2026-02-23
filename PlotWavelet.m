function PlotWavelet(W, cellIndex)
    figure;
    surface(1:size(W.cfs, 2), W.freqs, abs(W.cfs));
    shading interp
    axis tight
    set(gca, 'yscale', 'log');
    title(['Wavelet transform for cell ' num2str(cellIndex)]);
    colorbar
end

