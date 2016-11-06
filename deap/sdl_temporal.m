%%
% clear
% channel = 15;
% datFiles = dir('./s*.mat');
% % Fun, anger, sad
% % movies = [1, 3, 11, 12, 14,...
% %     16, 22, 23, 24, 25,...
% %     35, 37, 38, 39, 40];
% movies = [1, 3, 11, ...
%     16, 22, 23,...
%     35, 37, 38];
% vecPerClass = 3;
% userCount = size(datFiles, 1);
% X = zeros(8064, size(movies, 2)*userCount);
% for fi = 1:userCount
%     disp(['User file:', datFiles(fi).name])
%     load(datFiles(fi).name);
%     col = 0;
%     for m = movies
%         loc = mod(col, vecPerClass) + 1 + vecPerClass*(fi - 1) + floor(col/vecPerClass)*userCount*vecPerClass;
%         X(:, loc) = data(m, channel, :);
%         col = col+1;
%     end
%     clear data labels
% end
% X = X - repmat(mean(X,2), 1, size(X,2));

clear
X = load('loadedSmallOrder.mat');
X = X.X;
pwr = @(s) sum(s.*s, 1);
mag = @(v) sum(v.^2, 1);
%% Plot Data
plot(X(1:100,1:5))
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

[D, W] = ksvd(X, 200, 20);
%%
figure
clf
imagesc(W~=0)
%% Presence of emotion information in dictionary atoms
happy = [];
sad = [];
exciting = [];
users=  32;
for i = 1:users
    for col = 1:3
        happy = [happy;9*(i - 1) + col];
    end
end

for i = 1:users
    for col = 4:6
        sad = [sad;9*(i - 1) + col];
    end
end

for i = 1:users
    for col = 7:9
        exciting = [exciting;9*(i - 1) + col];
    end
end

clf
atomPresence = (W(:, happy) > 0);
atomPresence = sum(atomPresence, 2);
subplot(3, 1, 1)
plot(atomPresence)
find(atomPresence == max(atomPresence))
title('Frequency of Dictionary atoms for signals(Happy feeling)')

atomPresence = (W(:, sad) > 0);
atomPresence = sum(atomPresence, 2);
subplot(3, 1, 2)
plot(atomPresence)
find(atomPresence == max(atomPresence))
title('Frequency of Dictionary atoms for signals(Sad)')

atomPresence = (W(:, exciting) > 0);
atomPresence = sum(atomPresence, 2);
subplot(3, 1, 3)
plot(atomPresence)
find(atomPresence == max(atomPresence))
title('Frequency of Dictionary atoms for signals(Excited)')

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

%% Compare Time Signals

X_a = D*W;
clf
c = 65;
plot(X(:,c), 'r')
hold on
plot(X_a(:, c), 'b')

%% Plot dictionary atom
plot(D(:,250))
%%
figure()
m = mean(X')';
plot(m)