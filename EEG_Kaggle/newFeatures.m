f1 = load('data/train_1/1_10_0.mat');
%% 
dataS = f1.dataStruct.data;
part1 = dataS(1:4096, 1);

%% Spectral Entropy
xf = abs(fft(part1));
figure(2)
plot(xf)
xf = xf.^2;
xf = xf/sum(xf);
ent = -sum(xf.*log2(xf));

%% Spectrum correlation matrix across channels
close all
xf = zeros(4096, 16);
for c = 1:16
    xf(:,c) = abs(fft(dataS(1:4096,c)));
end
spec_corrmat = corrcoef(xf);
figure
imagesc(spec_corrmat, [-1, 1]);
%% 16x16 correlation matrix of time series across channels
close all
dataC = num2cell(dataS, 1);
figure
corrmat = corrcoef(dataS);
imagesc(corrmat)