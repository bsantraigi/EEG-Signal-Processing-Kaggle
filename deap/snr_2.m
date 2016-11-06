function [ SNR ] = snr_2( actual, approx )
%SNR_2 Summary of this function goes here
%   Detailed explanation goes here
pwr = @(s) sum(s.*s, 1);
n = approx - actual;

SNR = pwr(approx)./pwr(n);

SNR = (10*log10(SNR));

end

