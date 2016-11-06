%% Load Matrix Format 1
clear
close all
channel = 15;
datFiles = dir('./s*.mat');
% Fun, anger, sad
% movies = [1, 3, 11, 12, 14,...
%     16, 22, 23, 24, 25,...
%     35, 37, 38, 39, 40];
movies = [1, 3,...
    16, 22,...
    35, 37];
vecPerClass = 2;
% userCount = size(datFiles, 1);
userCount = 6;
n_samples = 8000;
X = zeros(n_samples, size(movies, 2)*userCount);
for fi = 1:userCount
    disp(['User file:', datFiles(fi).name])
    load(datFiles(fi).name);
    col = 0;
    for m = movies
        p = mod(col, vecPerClass) + 1 + vecPerClass*(fi - 1) + floor(col/vecPerClass)*userCount*vecPerClass;
        locS = (p-1)*32 + 1;
        locE = locS + 32 - 1;
        fprintf('%d, %d\n', locS, locE)
        X(:, locS:locE) = squeeze(data(m, 1:32, 1:n_samples))';
        col = col+1;
    end
    clear data labels
end
X = X - repmat(mean(X,2), 1, size(X,2));

%% Load Matrix Format 2
clear
close all
datFiles = dir('./s*.mat');
movies = 1;
userCount = 3;
n_samples = 1000;
X = zeros(32, userCount*n_samples);
for fi = 1:userCount
    disp(['User file:', datFiles(fi).name])
    load(datFiles(fi).name);
    s = ((fi-1)*n_samples + 1);
    X(:, s:(s + n_samples - 1)) = squeeze(data(movies, 1:32, 1:n_samples));
    clear data labels
end
X = X - repmat(mean(X,2), 1, size(X,2));

%% Get energy mat from X
pwr2 = @(s) sqrt(sum(s.*s, 2));
packetSize = 100;
E = zeros(32, userCount*n_samples/packetSize);
for i = 1: size(E, 2)
    f = (i-1)*packetSize + 1;
    E(:, i) = pwr2(X(:, f:(f+packetSize - 1)));
end
%% 
% clear
% X = load('loadedSmallOrder.mat');
% X = X.X;
pwr = @(s) sum(s.*s, 1);
mag = @(v) sum(v.^2, 1);
%% Plot Data
plot(X(:,1:100))
%% Optional
% X = fftshift(fft(X, 8192));
% X = 10*log10((X.*conj(X)).^0.5);
%% Decide possible rank of X
clf

[u,s,v] = svd(X);
s = max(s);
fprintf('Singular value range of X: %f, %f\n', min(s), max(s));
plot(s)

%% Run KSVD

% [D, W] = ksvd(X, 100, 20, 4);
[D, W] = ksvd(E, 100, 20, 4);
%%
figure
clf
imagesc(W~=0)

%% Correlation b/w the Dictionary atoms
corrmat = zeros(size(D, 2), size(D, 2));
for i = 1:size(D, 2)
    disp(i)
    for j = i:size(D, 2)
        corrmat(i, j) = xcorr(...
            (D(:,i) - mean(D(:,i)))/std(D(:,i)),...
            (D(:,j) - mean(D(:,j)))/std(D(:,j)), 0);
        corrmat(j, i) = corrmat(i, j);
    end
end
corrmat = corrmat/size(D,1);
%% Find basis vectors with high correlations
[is, js] = find(abs(corrmat) > 0.8);
rcc = [is, js];
rcc = rcc(is ~= js, :);
ids = unique(rcc(:))';

eqvMat = zeros(size(D, 2), size(D, 2));

for a = rcc'
    eqvMat(a(1), a(2)) = 1;
    eqvMat(a(2), a(1)) = 1;
end

while sum(sum(eqvMat, 1) > 1) > 0
    for a = rcc'
        eqvMat(a(1), eqvMat(a(2), :) > 0) = 1;
        eqvMat(a(2), :) = 0;
    end
end

toRemove = 0;
for i = 1:size(D, 2)
    if(sum(eqvMat(i, :)) > 0)
        u = unique([i, find(eqvMat(i,:) > 0)])';
        toRemove = toRemove + size(u, 1) - 1;
    end
end
fprintf('You can remove %d atoms.\n', toRemove)
% clear eqvMat rcc is js a ids

%%

%%
% 22 dB is good for freq. domain
X_a = D*W;
n = X_a - X;

SNR = pwr(X)./pwr(n);

SNR = (10*log10(SNR));
disp(sprintf('Median: %0.2f Max: %0.2f Min: %0.2f', median(SNR), max(SNR), min(SNR)))
c = find(SNR == min(SNR));
%% Compare Time Signals
figure()
X_a = D*W;
clf
plot(X(:,c), 'r')
hold on
plot(X_a(:, c), 'b')

%% Plot dictionary atom
clf
l = size(D, 2);
l = 25;
g = ceil(sqrt(l));
for i = 1:l
    subplot(g, g, i)
    plot(D(:,i))
end
%%
figure()
m = mean(X')';
plot(m)