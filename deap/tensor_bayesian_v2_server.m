clearvars
close all
%%
load('s01.mat')
p1_data = data;
p1_labels = labels;

clear data labels

load('s10.mat')
p1_data2 = data;
p1_labels2 = labels;

clear data labels
%% Draw fft of data
Fs = 128;
NFFT = 1024; % Next power of 2 from length of y
L = 1024;
Y = squeeze(fft(p1_data(2,1,1:1024),NFFT)/L);
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(1)
clf
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
%%
close all
C = 32; % First dimension of Y
% Taking 1 sec of data, Fs = 128 Hz, from deap site
T = 68; % Second dimension of Y
N = 10; % Third dimension/depth of Y
Y = zeros(C,T,N);
%%
clc
n_movies = N;
lastM = 0;
t = 0;
for p = 1:N
    m = ceil(p/(N/n_movies))
    if lastM == m
        % increment in time axis
        t = t + T;
        Y(:,:,p) = p1_data2(m, 1:C, (t + 1):(t + T));        
    else
        t = 0;
        Y(:,:,p) = p1_data2(m, 1:C, 1:T);
    end
    lastM = m;
end
%%
Y = MultiUserData([1, 2], C, T, N);
%%
figure(1)
clf
for i = 1:N
    imagesc(Y(:,:,i));
    pause
end
%%
K = 40; % Third dimension/depth of Y
D = zeros(C, T, K);
S = zeros(N, K);
gm_d = zeros(C, 1);
gm_s = zeros(N, K);
gm_n = zeros(C, N);

%% FIRST SAMPLES
clc
% Initialize gm_d
alpha_d = 1.6;
beta_d = 0.07;
for i = 1:C
    gm_d(i) = gamrnd(alpha_d, 1/beta_d);
end

% Initialize gm_s
alpha_s = 1.5;
beta_s = 0.06;
for i = 1:N
    for j = 1:K
        gm_s(i, j) = gamrnd(alpha_s, 1/beta_s);        
    end
end


% Initialize gm_n
alpha_n = 3;
beta_n = 0.3;

for i = 1:C
    for j = 1:N
        gm_n(i, j) = gamrnd(alpha_n, 1/beta_n);
    end
end

constants = struct();
constants.alpha_d = alpha_d;
constants.beta_d = beta_d;
constants.alpha_s = alpha_s;
constants.beta_s = beta_s;
constants.alpha_n = alpha_n;
constants.beta_n = beta_n;

% Initialize D
for j = 1:T
    for i = 1:C
        for k = 1:K
            if j == 1
                D(i, j, k) = normrnd(0, sqrt(1/gm_d(i)));
            else
                D(i, j, k) = normrnd(D(i, j - 1, k), sqrt(1/gm_d(i)));
            end
        end
    end
end

% Initialize S
for p = 1:N
    for k = 1:K
        S(p, k) = normrnd(0, 1/sqrt(gm_s(p, k)));
    end
end
%%
figure(1)
clf
imagesc(D(:,:,1));
%% START SAMPLING
iters = 300;
burn = 150;
meanD = zeros(C,T,K);
meanS = zeros(N,K);
for it = 1:iters
    disp(['Iteration:' num2str(it)])
    if it > burn
        meanD = meanD + D;
        meanS = meanS + S;
    end
    tic
    [D, S, gm_d, gm_s, gm_n] = GibbsSampleNextTensor_V2(Y, D, S, gm_d, gm_s, gm_n, constants, true);
    toc
end
%%
meanD = meanD/(iters - burn);
meanS = meanS/(iters - burn);
Y_approx = zeros(C, T, N);
for j = 1:T
    Y_approx(:,j,:) = squeeze(meanD(:,j,:))*meanS';
end
% save('fullWSpace_n40');
display('Mean Calculated...')
%%
close all
figure(1)
for c = 1:32
    clf
    hold on
    plot(Y(c, :, 1), 'r')
    plot(Y_approx(c, :, 1), 'b')
    hold off
    legend('Actual', 'Approx');
    pause
end
%% Check if approximation is good
close all
figure(1)
clf

for k = 1:N
    clf
    subplot(2, 1, 1)
    M = squeeze(Y(:,:,k));
    imagesc(M), colorbar
    xlabel('Time n')
    ylabel('Channel i(1~32)')
    title('Actual Data')
    m1 = min(M(:));
    m2 = max(M(:));
    subplot(2, 1, 2)
    imagesc(Y_approx(:,:,k)), colorbar
    xlabel('Time n')
    ylabel('Channel i(1~32)')
    title('Reconstructed Data')
    pause
end

%%
figure(1)
clf
for k = 1:K
    imagesc(meanD(:,:,k)), colorbar
    pause
end
%%
figure(1)
clf
imagesc(meanS'), colorbar

%%
figure(1)
clf
for k = 1:N
    imagesc(Y(:,:,k)), colorbar
    pause
end

