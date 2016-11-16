clear
clc
for user = 2:2
%     user = 1;

    % Training to be done using all available samples
    [safeTable_0, fselect0] = getSafeList(user, 0);
    [safeTable_1, fselect1] = getSafeList(user, 1);
    p0 = size(safeTable_0, 1)/(size(safeTable_0, 1) + size(safeTable_1, 1));
    p1 = size(safeTable_1, 1)/(size(safeTable_0, 1) + size(safeTable_1, 1));
    %%
    maxList_0 = zeros(256,1);
    minList_0 = 9e9*ones(256,1);
    maxList_1 = zeros(256,1);
    minList_1 = 9e9*ones(256,1);

    fcount = 5;
    gap = 1024;
    FS0 = zeros(length(1:gap:(240000 - 8192))*fcount, 256);
    FS1 = zeros(length(1:gap:(240000 - 8192))*fcount, 256);
    r = 1;

    for fi = 1:fcount
        fprintf('Learning from File #%d\n', fi);
        f0 = load(safeTable_0.path{fi});
        x0 = f0.dataStruct.data;
        f1 = load(safeTable_1.path{fi});
        x1 = f1.dataStruct.data;

        window = (zeros(240000, 1) == 1);
        window(1:8192) = 1;

        for shift = 1:gap:(240000 - 8192)    
    %         fprintf('Shift: %d\n', shift);
            window = circshift(window, gap);
            fs0 = fExtractor_v2(x0(window, :));
            fs1 = fExtractor_v2(x1(window, :));
            fs0 = fs0(:);
            fs1 = fs1(:);

            FS0(r, :) = fs0(:);
            FS1(r, :) = fs1(:);
            r = r+1;

            b_0 = minList_0 > fs0;
            minList_0(b_0) = fs0(b_0);
            b_0 = maxList_0 < fs0;
            maxList_0(b_0) = fs0(b_0);

            b_1 = minList_1 > fs1;
            minList_1(b_1) = fs1(b_1);
            b_1 = maxList_1 < fs1;
            maxList_1(b_1) = fs0(b_1);
        end
    end
    %%
    binSize = 50;
    histMat_0 = zeros(256, binSize);
    histMat_1 = zeros(256, binSize);
    histEdges_0 = zeros(256, binSize + 1);
    histEdges_1 = zeros(256, binSize + 1);

    % bins_0 = zeros(256, 21);
    % bins_1 = zeros(256, 21);
    % for f = 1:256
    %     bins_0(f, :) = minList_0(f):(maxList_0(f) - minList_0(f))/20:maxList_0(f);
    %     bins_1(f, :) = minList_1(f):(maxList_1(f) - minList_1(f))/20:maxList_1(f);
    % end

    % FS0 = zeros(length(1:gap:(240000 - 8192)), 256);
    % FS1 = zeros(length(1:gap:(240000 - 8192)), 256);
    % r = 1;
    % 
    % window = (zeros(240000, 1) == 1);
    % window(1:8192) = 1;
    % for shift = 1:gap:(240000 - 8192)
    %     fprintf('Shift: %d\n', shift);
    %     window = circshift(window, gap);
    %     fs0 = fExtractor_v2(x0(window, :));
    %     fs1 = fExtractor_v2(x1(window, :));
    %     FS0(r, :) = fs0(:);
    %     FS1(r, :) = fs1(:);
    %     r = r+1;
    % end

    % FIXME(TODO): REMOVE OUTLIERS - USE M.A.D
    madThresh = 3;
    for f = 1:256
        fs0 = FS0(:, f);
        m0 = median(fs0);
        mad0 = median(abs(fs0 - m0));
        fs0(abs((fs0 - m0)/mad0) > madThresh)= [];
        
        fs1 = FS1(:, f);
        m1 = median(fs1);
        mad1 = median(abs(fs1 - m1));
        fs1(abs((fs1 - m1)/mad1) > madThresh)= [];
        
        [histMat_0(f, :), histEdges_0(f, :)] = histcounts(fs0, binSize);
        histMat_0(f, :) = histMat_0(f, :)/sum(histMat_0(f, :));
        [histMat_1(f, :), histEdges_1(f, :)] = histcounts(fs1, binSize);
        histMat_1(f, :) = histMat_1(f, :)/sum(histMat_1(f, :));
    end    
    
    figure(2)
    clf
    pic = 0;
    for f = 1:100
        pic = pic + 1;
        subplot(10, 10, pic);
        hold on
        plot(histEdges_0(f, 1:end-1),histMat_0(f, :))
        plot(histEdges_1(f, 1:end-1), histMat_1(f, :))
        hold off
    end

    %%
    fselect = GetTestFiles(user);
    resultMap = containers.Map('KeyType','char','ValueType','double');
    PF = 0.01;
    PD = 0.98;
    tau_0 = log((1-PD)/(1-PF));
    tau_1 = log(PD/PF);

    for fi = 2:length(fselect)
        fname = fselect{fi};
        ftest = load(fname);
        xtest = ftest.dataStruct.data;

        window = (zeros(240000, 1) == 1);
        window(1:8192) = 1;
        L = 0;
        logL0 = 0;
        logL1 = 0;
        for shift = 1:100:(240000 - 8192)
    %         fprintf('Shift: %d\n', shift);
            window = circshift(window, 100);
            fs = fExtractor_v2(xtest(window, :));
            fs = fs(:);
            [loglr, like0, like1] = CalcLogLR(fs, histMat_0, histEdges_0, histMat_1, histEdges_1);
            logL0 = logL0 + like0;
            logL1 = logL1 + like1;
            L = L + loglr;

            if(L <= tau_0)
                fprintf('%s: H0\n', fname(17:end));
                resultMap(fname(17:end)) = p1/(p1 + p0*exp(logL0 - logL1));
                fprintf('LRT: %f, P(H1/D): %f\n', L, p1/(p1 + p0*exp(logL0 - logL1)));
                break;
            end
            if(L >= tau_1)
                fprintf('%s: H1\n', fname(17:end));
                resultMap(fname(17:end)) = p1/(p1 + p0*exp(logL0 - logL1));
                fprintf('LRT: %f, P(H1/D): %f\n', L, p1/(p1 + p0*exp(logL0 - logL1)));
                break;
            end        
        end
    end

    %% Import 0.54 AUC
    filename = 'seqSol.csv';
    delimiter = ',';
    startRow = 2;
    formatSpec = '%s%f%[^\n\r]';
    fileID = fopen(filename,'r');

    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

    fclose(fileID);

    rsol_File = dataArray{:, 1};
    rsol_Class = dataArray{:, 2};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    filename = 'MatFiles_7Bands_new\MetaFile_Test.csv';
    delimiter = ',';
    formatSpec = '%*s%s%f%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    fNames = dataArray{:, 1};
    fIndex = dataArray{:, 2};
    clearvars filename delimiter formatSpec fileID dataArray ans;
    %% Import meta for user 2
    % pData = cell(size(yfit, 1)/29,2);
    j = 1;
    for k = 1:length(rsol_File)
        if ~isKey(resultMap, rsol_File(k))
            resultMap(rsol_File{k}) = rsol_Class(k);
        end
    end
    % needle = sprintf('%d_', user)
    % for k = 1:length(fNames)
    %     fName = strtrim(fNames{k});
    %     if(strfind(fName, needle) > 0)
    %         i = fIndex(k);
    % %         fprintf('Replace "%s": %f with %f\n', fName, resultMap(fName), yfit(i));
    %         resultMap(fName) = yfit(i);
    % %         pData{j, 1} = fName;
    % % y        pData{j, 2} = yfit(i);
    %         j = j+1;
    %     end
    % end
    solFile = fopen('seqSol.csv', 'w');
    fprintf(solFile, 'File,Class\n');
    for i = 1:length(fNames)
        fName = strtrim(fNames{i});
        fprintf(solFile, '%s,%f\n', fName, resultMap(fName));
    end
    fclose(solFile);
end












