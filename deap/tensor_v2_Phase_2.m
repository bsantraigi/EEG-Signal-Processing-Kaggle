%% PART 2 - Estimate Selection matrix for Test Data
% Select data from user 1 and user 2
% Select data other the training data
% Tensor D is fixed in phase 2

%% Create Data Matrix, Y
clear all
clearvars
load('fullWSpace_MultiUser.mat')

outFileName_S = sprintf('Estimate_SMat_MultiUser.mat')
%%
close all
C = 32; % First dimension of Y
% Taking 1 sec of data, Fs = 128 Hz, from deap site
T = 120; % Second dimension of Y
N = 40; % Third dimension/depth of Y
Y = MultiUserData([1, 3, 5], C, T, N);

%% Create variables - Set 2

% K = ~, Load from memory
D = zeros(C, T, K);
S = zeros(N, K);
gm_d = zeros(C, 1);
gm_s = zeros(N, K);
gm_n = zeros(C, N);

%% FIRST SAMPLES
% Initialize gm_d
alpha_d = 1.6;
beta_d = 0.07;
for i = 1:C
    gm_d(i) = gamrnd(alpha_d, 1/beta_d);
end

% Initialize gm_s
alpha_s = 300;
beta_s = 300;
for i = 1:N
    for j = 1:K
        gm_s(i, j) = gamrnd(alpha_s, 1/beta_s);        
    end
end


% Initialize gm_n
alpha_n = 300;
beta_n = 300;

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

% D IS FIXED IN PHASE 2
D = meanD;

% Initialize S
for p = 1:N
    for k = 1:K
        S(p, k) = normrnd(0, 1/sqrt(gm_s(p, k)));
    end
end

%% START SAMPLING
iters = 400;
burn = 200;
meanS = zeros(N,K);
for it = 1:iters
    disp(['Iteration:' num2str(it)])
    if it > burn
        meanS = meanS + S;
    end
    tic
    [~, S, gm_d, gm_s, gm_n] = GibbsSampleNextTensor_V2(Y, D, S, gm_d, gm_s, gm_n, constants, false);
    toc
end

meanS = meanS/(iters - burn);
save(outFileName_S, 'meanS');


%% Check if approximation is good
Y_approx = zeros(C, T, N);
for j = 1:T
    Y_approx(:,j,:) = squeeze(meanD(:,j,:))*meanS';
end

figure(1)
clf
imagesc(meanS'), colorbar
pause

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










