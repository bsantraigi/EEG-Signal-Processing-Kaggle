%% Load Data
clear all
user = 3;
dinter = load(sprintf('MatFiles_7Bands_&Meta/train_%d[interictal].mat', user));
dpre = load(sprintf('MatFiles_7Bands_&Meta/train_%d[preictal].mat', user));
li = size(dinter.newData, 2);
lp = size(dpre.newData, 2);
prob_preictal = lp/(lp + li);
li = min(li, lp*4);
data = [dinter.newData(:, 1:li), dpre.newData(:, 1:lp)]';

labels = [zeros(li, 1); ones(lp, 1)];
dataCombo = [data labels];

% labels = [ones(li, 1); -ones(lp, 1)];
% labels = 1.5*[labels, -labels];
labels = [zeros(li, 1); ones(lp, 1)];
%% Create ANN
fAlgo = sprintf('NN_trained_user_%d.mat', user)
if exist(fAlgo) == 2
    load(fAlgo);
else
    net = fitnet([20 20 20]);
    net = train(net, data', labels', 'useParallel','yes');
    save(sprintf('NN_trained_user_%d', user), 'net', '-v7.3');
end
%% Plot Test Data - Test out trained ANN
figure(1)
labP = net(data');
perform(net, labP, labels')
% labP(2,:) = labP(2,:)*prob_preictal;
% labP(1,:) = labP(1,:)*(1-prob_preictal);
% labP = softmax(labP)';

for i = 1:size(labP, 1)/29
%     labP(((i - 1)*29 + 1):(i*29), 2) = median(labP(((i - 1)*29 + 1):(i*29), 2));
    labP(((i - 1)*29 + 1):(i*29), 1) = median(labP(((i - 1)*29 + 1):(i*29), 1));
end

labP(labP < 0) = 0;
labP(labP > 1) = 1;

figure(1)
clf
subplot(2, 1, 1)
hold on
plot(labels(:,1))
plot(labP)
axis([0 2*li -1.5 1.5])
title('Interictal')
hold off

% subplot(2, 1, 2)
% hold on
% plot(labels(:,2)/3 + 0.5)
% plot(labP(:,2))
% axis([0 2*li -1.5 1.5])
% title('Preictal')
% hold off

%% Prediction Step
sprintf('MatFiles_7Bands_&Meta/test_%d.mat', user)
dtest = load(sprintf('MatFiles_7Bands_&Meta/test_%d.mat', user));
dtest = [dtest.newData];
yPredict = net(dtest);
% yPredict(2,:) = yPredict(2,:)*prob_preictal;
% yPredict(1,:) = yPredict(1,:)*(1-prob_preictal);
% py_predict = softmax(yPredict)';

py_predict = yPredict;
py_predict(py_predict < 0) = 0;
py_predict(py_predict > 1) = 1;

for i = 1:size(py_predict, 1)/29
%     py_predict(((i - 1)*29 + 1):(i*29), 2) = median(py_predict(((i - 1)*29 + 1):(i*29), 2));
    py_predict(((i - 1)*29 + 1):(i*29), 1) = median(py_predict(((i - 1)*29 + 1):(i*29), 1));
end

figure(1)
clf

subplot(2, 1, 1)
plot(py_predict)
title('Interictal')

% subplot(2, 1, 2)
% plot(py_predict(:,2))
% title('Preictal')
% hold off

yfit = py_predict;

%% Import 0.54 AUC
filename = '/home/bishal/Documents/MTP/EEG_Kaggle/NNSol.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

rsol_File = dataArray{:, 1};
rsol_Class = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
filename = '/home/bishal/Documents/MTP/EEG_Kaggle/MatFiles_7Bands_&Meta/MetaFile_Test.csv';
delimiter = ',';
formatSpec = '%*s%s%f%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
fNames = dataArray{:, 1};
fIndex = dataArray{:, 2};
clearvars filename delimiter formatSpec fileID dataArray ans;
%% Import meta for user
% pData = cell(size(yfit, 1)/29,2);
resultMap = containers.Map('KeyType','char','ValueType','double');
j = 1;
for k = 1:length(rsol_File)
    resultMap(rsol_File{k}) = rsol_Class(k);
end
needle = sprintf('%d_', user)
for k = 1:length(fNames)
    fName = strtrim(fNames{k});
    if(strfind(fName, needle) > 0)
        i = fIndex(k);
%         fprintf('Replace "%s": %f with %f\n', fName, resultMap(fName), yfit(i));
        resultMap(fName) = yfit(i);
%         pData{j, 1} = fName;
% y        pData{j, 2} = yfit(i);
        j = j+1;
    end
end
solFile = fopen('NNSol.csv', 'w');
fprintf(solFile, 'File,Class\n');
for i = 1:length(fNames)
    fName = strtrim(fNames{i});
    fprintf(solFile, '%s,%0.12f\n', fName, resultMap(fName));
end
fclose(solFile);