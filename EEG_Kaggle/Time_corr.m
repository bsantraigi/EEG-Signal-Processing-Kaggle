user = 2;
f1 = load(sprintf('data/train_%d/%d_12_0.mat', user, user));
f1 = f1.dataStruct;
f2 = load(sprintf('data/train_%d/%d_10_1.mat', user, user));
f2 = f2.dataStruct;
%%
train_f1 = load('C:\Users\Reincarnated\Documents\MTP\EEG_Kaggle\MatFiles_7Bands_&Meta\train_3[interictal].mat');
train_f2 = load('C:\Users\Reincarnated\Documents\MTP\EEG_Kaggle\MatFiles_7Bands_&Meta\train_3[preictal].mat');