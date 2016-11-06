function [newData] = preprocessTrainData(user, type2, fileLimit)
    n = 1; % Downsample Factor
    Fs = 400/n;
    
    % Standard Frequency bands - alpha, beta, delta and gamma
    bands = [
        [0, 4, 7, 12, 40, 70, 110];
        [4, 7, 12, 40, 70, 110, 200]
    ];
    
    % N-DFT: value of 'N'
    pt = 8192;

    bands = floor(pt/Fs*bands);
    bands(1,:) = bands(1,:)+1;
    outerColIndex = 1;
    
    % Data ordering
    % -- user 1, 1 hour, interictal ||| user 1, 1 hour, preictal ----
    if strcmp(type2, 'interictal')
        type = 0;
    elseif strcmp(type2, 'preictal')
        type = 1;
    end
%     needle = sprintf('_%d.mat', type);
%     fold = sprintf('train_%d/', user);
    fselect = {};
    
    h = waitbar(0, '', 'Name', 'Looking for files...');
    steps = fileLimit;
    
    % Enumerate through all possible sequence ids 
    % Choose the ones that actually exists
    for i = 1:fileLimit
        fname = sprintf('data/train_%d/%d_%d_%d.mat', user, user, i, type);
        waitbar(i / steps, h, sprintf('data/train-%d/%d-%d-%d.mat', user, user, i, type));
        if(exist(fname) == 2)      
            fselect = [fselect {fname}];
        end
    end
    
    % P.S. - Each data file (containing 10 min. data) is converted to
    % 29 seperate feature columns which are uniformly distributed
    % in time i.e. [0 ~ 10/29, 10/29 ~ 20/29, ..., 280/29 ~ 10] mins.
    newData = zeros(80, 29*length(fselect));
    
    close(h);    
    h = waitbar(0, '', 'Name', 'Generating Data Matrix...');
    steps = length(fselect);
    
    % Create a meta file with the details about the location 
    % of the files in the combined data matrix
    metaFileID = fopen('MetaFile_Train.csv', 'a');
    
    % Now load data files one by one
    % Create 29 feature columns for each of them
    % and then append them to the newData matrix
    
    for i = 1:length(fselect)
        fname = fselect{i};
        f1 = load(fname);
        for ch = 1:16
            % Pickup channel data
            x_interictal = f1.dataStruct.data(:, ch);
            xd_interictal = downsample(x_interictal, n);
            colIndex = outerColIndex;
            
            % Group time series into 8192 point
            % segments - Calculate fft of all the segments
            % Calculate band energies
            for t = 1:pt:length(xd_interictal)
                if(t+pt > length(xd_interictal))
                    break
                end
                ff = fft(xd_interictal(t:(t+pt-1)));
                ff = abs(ff)/1e5;

                for bx = 1:size(bands, 2)
                    newData((ch-1)*size(bands, 2) + bx, colIndex) = norm(ff(bands(1, bx):bands(2, bx)));
                end
                % Append new features from index (size(bands, 2)*16ch + 1)
                
                % Increment column for next time segment
                colIndex = colIndex + 1;
            end
        end
        
        colIndex = outerColIndex;        
        xd_interictal = f1.dataStruct.data;
        for t = 1:pt:size(xd_interictal, 1)
            if(t+pt-1 > size(xd_interictal, 1))
                break
            end            
            fnew = fExtractor(xd_interictal(t:(t+pt-1), :));
            fnew = fnew(:);
            % Append new features from index (size(bands, 2)*16ch + 1)
            enter = (size(bands, 2)*16 + 1);
            newData(enter:(enter + length(fnew) - 1), colIndex) = fnew;
            % Increment column for next time segment
            colIndex = colIndex + 1;
        end
        
        col_s = outerColIndex;
        col_f = colIndex - 1;
        fsplit = strsplit(fname, '/');
        fprintf(metaFileID, '%s, %s, %d, %d\n', fsplit{2}, fsplit{3}, col_s, col_f);
        
        outerColIndex = outerColIndex + 29;
        waitbar(i / steps, h, sprintf('data/train-%d/%d-%d-%d.mat', user, user, i, type));
    end
    fclose(metaFileID);
    close(h);

end