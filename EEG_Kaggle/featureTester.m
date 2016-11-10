%% Amplitude histogram
%% Measure min and max
for user = 1:3
    [~, fselect0] = getSafeList(user, 0);
    [~, fselect1] = getSafeList(user, 1);
    minamp = 10000;
    maxamp = -10000;
    for i = 1:length(fselect0)
        f = fselect0{i};
        f = load(f);
        minamp = min(min(min(f.dataStruct.data)), minamp);
        maxamp = max(max(max(f.dataStruct.data)), maxamp);
    end
    for i = 1:length(fselect1)
        f = fselect1{i};
        f = load(f);
        minamp = min(min(min(f.dataStruct.data)), minamp);
        maxamp = max(max(max(f.dataStruct.data)), maxamp);
    end
    fprintf('Range for user %d: (%d, %d)\n', user, minamp, maxamp);
%     Range for user 1: (-2.102788e+03, 1.991541e+03)
%     Range for user 2: (-1.637427e+03, 9.683459e+02)
%     Range for user 3: (-1.424802e+03, 8.241125e+02)
end
%% Histogram bins
left = 0;
right = 0;
k = 1.5;
edges = [0];
for i = 1:20    
    left = left - i^k;
    right = right + i^k;
    edges = [left edges right];
end
%% Load filelist
user = 2;
[~, fselect0] = getSafeList(user, 0);
[~, fselect1] = getSafeList(user, 1);
%% Load File
% f0 = load(sprintf('data/train_%d/%d_70_0', user, user));
% f1 = load(sprintf('data/train_%d/%d_2_1', user, user));
f0 = load(fselect0{15});
f1 = load(fselect1{50});
% Load channel
y0 = f0.dataStruct.data(:,1);
y1 = f1.dataStruct.data(:,1);

%% Histogram
figure(1)
clf
hold on
histogram(y1, edges)
histogram(y0, edges)
hold off

%% Get feature matrix with label column attached

left = 0;
right = 0;
k = 1.3;
edges = [0];
for i = 1:20
    left = left - i^k;
    right = right + i^k;
    edges = [left edges right];
end

for user = 1:3
    [~, fselect0] = getSafeList(user, 0);
    [~, fselect1] = getSafeList(user, 1);
    fd = length(edges);
    feat3DMat = zeros(length(fselect0) + length(fselect1),...
         fd*16 + 1);
    h = waitbar(0, '', 'Name', 'Extracting Juice...');
    steps = length(fselect0) + length(fselect1);
    for i = 1:length(fselect0)
        waitbar(i / steps, h, sprintf('Progress %0.2f%%', 100*i / steps));
        f = fselect0{i};
        f = load(f);
        for c = 1:16
            feat3DMat(i, ((c-1)*fd + 1):(c*fd) ) = histc(f.dataStruct.data(:, c), edges);
        end
        feat3DMat(i, end) = 0;
    end
    offset = length(fselect0);
    for i = 1:length(fselect1)
        waitbar((offset + i) / steps, h, sprintf('Progress %0.2f%%', 100*(offset + i) / steps));
        f = fselect1{i};
        f = load(f);
        for c = 1:16
            feat3DMat(offset + i, ((c-1)*fd + 1):(c*fd) ) = histc(f.dataStruct.data(:, c), edges);
        end
        feat3DMat(offset + i, end) = 1;
    end
    close(h);
    save(sprintf('feat3DMat_%d', user), 'feat3DMat', '-v7.3');
end

%%
rng(1)

SVMModel = fitcsvm(feat3DMat(:, 1:end-1),feat3DMat(:, end),'Standardize',true,'KernelFunction','Polynomial', 'PolynomialOrder', 2, ...
    'KernelScale','auto');
CVSVMModel = crossval(SVMModel);
classLoss = kfoldLoss(CVSVMModel)

%%
figure(3)
clf
subplot(2,1,1)
plot(y0)
subplot(2,1,2)
plot(y1)