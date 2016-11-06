%%
clear
chRanks = zeros(40,1);
for channel = 21:40
    datFiles = dir('./s*.mat');
    disp(['Channel:' num2str(channel)]) 
    % Fun, anger, sad
    % movies = [1, 3, 11, 12, 14,...
    %     16, 22, 23, 24, 25,...
    %     35, 37, 38, 39, 40];
    movies = [1, 3, 11, ...
        16, 22, 23,...
        35, 37, 38];
    vecPerClass = 3;
    userCount = 10;
    X = zeros(8064,vecPerClass*3*userCount);
    % userCount = size(datFiles, 1);
    
    for fi = 1:userCount
        disp(['User file:', datFiles(fi).name])
        load(datFiles(fi).name);
        col = 0;
        for m = movies
            loc = mod(col, vecPerClass) + 1 + vecPerClass*(fi - 1) + floor(col/vecPerClass)*userCount*vecPerClass;
            X(:, loc) = data(m, channel, :);
            col = col+1;
        end
        clear data labels
    end
    X = X - repmat(mean(X')', 1, size(X,2));

    % clear
    % X = load('loadedSmallOrder.mat');
    % X = X.X;
    % pwr = @(s) sum(s.*s, 1);
    % mag = @(v) sum(v.^2, 1);
    %% Optional
    % X = fftshift(fft(X, 8192));
    % X = 10*log10((X.*conj(X)).^0.5);
    %% Decide possible rank of X
    disp('Calculating SVD');
    [u,s,v] = svd(X);
    s = max(s);
    % fprintf('Singular value range of X: %f, %f\n', min(s), max(s));
    plot(s)
    chRanks(channel) = find(s>max(s)/20, 1, 'last');
%     clear u s v X
end