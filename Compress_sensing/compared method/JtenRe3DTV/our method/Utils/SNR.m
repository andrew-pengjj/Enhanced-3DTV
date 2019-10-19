function snr=SNR(sig0,sig1)
 noiseE = sum((sig0(:) - sig1(:)).^2);
%  signalE = sum((sig1(:)-mean(sig1(:)) ).^2);
 signalE = sum( (sig1(:)).^2);
 snr = 10*log10(signalE/noiseE);


end