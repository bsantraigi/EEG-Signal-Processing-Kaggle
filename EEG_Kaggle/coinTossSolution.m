sFileID = fopen('randSolution.csv', 'w');
fprintf(sFileID, 'File,Class\n');

fl = dir('data/test_1');
for i = 1:length(fl)
    if(strcmp(fl(i).name, '..') == 1 || strcmp(fl(i).name, '.') == 1)
        continue
    end
    fprintf(sFileID, '%s,%f\n', fl(i).name, rand);
end

fl = dir('data/test_2');
for i = 1:length(fl)
    if(strcmp(fl(i).name, '..') == 1 || strcmp(fl(i).name, '.') == 1)
        continue
    end
    fprintf(sFileID, '%s,%f\n', fl(i).name, rand);
end

fl = dir('data/test_3');
for i = 1:length(fl)
    if(strcmp(fl(i).name, '..') == 1 || strcmp(fl(i).name, '.') == 1)
        continue
    end
    fprintf(sFileID, '%s,%f\n', fl(i).name, rand);
end

fclose(sFileID);