f = fopen('MetaFile_Train.csv', 'w');
fclose(f);
newData = preprocessTrainData(1, 'preictal', 1500);
save('train_1[preictal]', 'newData', '-v7.3');
clearvars newData

newData = preprocessTrainData(2, 'preictal', 2500);
save('train_2[preictal]', 'newData', '-v7.3');
clearvars newData

newData = preprocessTrainData(3, 'preictal', 2500);
save('train_3[preictal]', 'newData', '-v7.3');
clearvars newData

newData = preprocessTrainData(1, 'interictal', 1500);
save('train_1[interictal]', 'newData', '-v7.3');
clearvars newData

newData = preprocessTrainData(2, 'interictal', 2500);
save('train_2[interictal]', 'newData', '-v7.3');
clearvars newData


newData = preprocessTrainData(3, 'interictal', 2500);
save('train_3[interictal]', 'newData', '-v7.3');
clearvars newData

%%
f = fopen('MetaFile_Test.csv', 'w');
fclose(f);
newData = preprocessTestData(1, 1800);
save('test_1', 'newData', '-v7.3');
clearvars newData

newData = preprocessTestData(2, 3000);
save('test_2', 'newData', '-v7.3');
clearvars newData

newData = preprocessTestData(3, 3000);
save('test_3', 'newData', '-v7.3');
clearvars newData













