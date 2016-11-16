function [ fs ] = fExtractor_v2( partialDat )
%FEXTRACTOR Summary of this function goes here
%   partialDat - 8192*16 dimensional data

% Shannon's Entropy at seven frequency bands
bands = [
        [0, 4, 7, 12, 40, 70, 110];
        [4, 7, 12, 40, 70, 110, 200]
    ];
Fs = 400;
pt = 8192;
bands = floor(pt/Fs*bands);
bands(1,:) = bands(1,:)+1;

fs = zeros(16, 16);

for c = 1:16
    fftMag = abs(fft(partialDat(:, c)));
    % Energy in each bands
    for bx = 1:size(bands, 2)         
         fs(c, 9+bx) = norm(fftMag(bands(1, bx):bands(2, bx)));
    end
    
    if sum(fftMag) ~= 0
        % Spectral Entropy
        fftMag = fftMag.^2;        
        fftMag = fftMag/sum(fftMag);
        fs(c, 1) = -sum(fftMag.*log2(fftMag));
        
        % Frequency Band Entropy
        probs = zeros(size(bands, 2), 1);
        for bx = 1:size(bands, 2)
             probs(bx) = norm(fftMag(bands(1, bx):bands(2, bx)))^2;
        end
        probs = probs/sum(probs);
        fs(c, 2) = -sum(probs.*log2(probs));
    end
end

% Spectrum correlation matrix across channels
% and eigen values(3)
xf = zeros(8192, 16);
for c = 1:16
    xf(:,c) = abs(fft(partialDat(:,c)));
end

if sum(sum(xf)) ~= 0
    spec_corrmat = corrcoef(xf);
    fs(:, 3) = eig(spec_corrmat);
end

% 16x16 correlation matrix of time series across channels
% and the eigen values(4) of time-corr matrix
if sum(sum(partialDat)) ~= 0
    time_corrmat = corrcoef(partialDat);
    fs(:, 4) = eig(time_corrmat);
end

% Skewness(5) and kurtosis(6) & Hjorth parameters(7,8 and 9)
for c = 1:16
    y = partialDat(:, c);
    m = mean(y);
    s = std(y);
    if s ~= 0
        fs(c, 5) = 1/8192*sum((y - m).^3)/s^3; % skewness
        fs(c, 6) = 1/8192*sum((y - m).^4)/s^4; % kertosis
    end
    [fs(c, 7), fs(c, 8), fs(c, 9)] = HjorthParameters(y); % HjorthParameters (7,8,9)
end



% dyadic frequency bands


end

