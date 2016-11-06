%% Load Metadata
filename = 'EEG_Kaggle/MetaFile_Train.csv';
delimiter = ',';
formatSpec = '%*s%s%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
TrainMeta = table(dataArray{1:end-1}, 'VariableNames', {'fname','start','finish'});
clearvars filename delimiter formatSpec fileID dataArray ans;

TrainMeta.index = TrainMeta.finish/29;
TrainMeta.user = cellfun(@(x) str2num(x(1)), TrainMeta.fname, 'UniformOutput', false);
TrainMeta.preictal = cellfun(@(x) x((end-4):end), TrainMeta.fname, 'UniformOutput', false);
TrainMeta.preictal = cellfun(@(x) str2num(x(1)), TrainMeta.preictal, 'UniformOutput', false);
%% Prepare Data for user
user = 3;
UserTable = TrainMeta(cellfun(@(x) x == user, TrainMeta.user), :);
fs = sprintf('EEG_Kaggle/user%d_Dictionaries_set1.mat', user);
fuser = load(fs);
train_data = fuser.S.*fuser.B;
train_data2 = fuser.S2.*fuser.B2;
train_labels = cell2mat(UserTable.preictal);
train_data = [train_data' train_data2' train_labels];
train_data = train_data(1:800, :);
% mix = randperm(800, 800);
train_data = train_data(mix, :);

%%




