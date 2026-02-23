function W = ComputeWavelet(y)
[cfs, freqs] = cwt(y, 'amor');
W.cfs = cfs;
W.freqs = freqs;
end