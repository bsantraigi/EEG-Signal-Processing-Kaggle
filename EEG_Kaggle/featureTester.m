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

left = 0;
right = 0;
k = 1.5;
edges = [0];
for i = 1:15
    left = left - i^k;
    right = right + i^k;
    edges = [left edges right];
end

figure(1)
clf
hold on
histogram(y1, edges)
histogram(y0, edges)
hold off
legend('Preictal', 'Interictal')

%%
figure(3)
clf
subplot(2,1,1)
plot(y0)
subplot(2,1,2)
plot(y1)