function [newData] = preprocessTestData(user, fileLimit)
    n = 1; % Downsample Factor
    Fs = 400/n;
    bands = [
        [0, 4, 7, 12, 40, 70, 110];
        [4, 7, 12, 40, 70, 110, 200]
    ];
    
    pt = 8192;

    bands = floor(pt/Fs*bands);
    bands(1,:) = bands(1,:)+1;
    outerColIndex = 1;
    % Data ordering
    % -- user 1, 1 hour, interictal ||| user 1, 1 hour, preictal ----

    fselect = {};
    
    h = waitbar(0, '', 'Name', 'Looking for files...');
    steps = fileLimit;
    
    for i = 1:fileLimit
        fname = sprintf('data/test_%d_new/new_%d_%d.mat', user, user, i);
        waitbar(i / steps, h, sprintf('data/test-%d-new/new-%d-%d.mat', user, user, i));
        if(exist(fname) == 2)
            fselect = [fselect {fname}];
        end
    end
    newData = zeros(80, 29*length(fselect));
    
    close(h);
    h = waitbar(0, '', 'Name', 'Generating Data Matrix...');
    steps = length(fselect);
    metaFileID = fopen('MetaFile_Test.csv', 'a');
    for i = 1:length(fselect)
        fname = fselect{i};
        f1 = load(fname);
        for ch = 1:16
            x_interictal = f1.dataStruct.data(:, ch);
            xd_interictal = downsample(x_interictal, n);
            colIndex = outerColIndex;
            for t = 1:pt:length(xd_interictal)
                if(t+pt > length(xd_interictal))
                    break
                end
                ff = fft(xd_interictal(t:t+pt));
                ff = abs(ff)/1e5;

                for bx = 1:size(bands, 2)
                    newData((ch-1)*size(bands, 2) + bx, colIndex) = norm(ff(bands(1, bx):bands(2, bx)));
                end
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
        
        waitbar(i / steps, h, sprintf('data/test-%d-new/new-%d-%d.mat', user, user, i));
    end
    fclose(metaFileID);
    close(h);

end