clear
close all
%%
load('s01.mat')
p1_data = data;
p1_labels = labels;

clear data labels

load('s02.mat')
p1_data2 = data;
p1_labels2 = labels;

clear data labels
%%
clc
close all
C = 32; % First dimension of Y
T = 20; % Second dimension of Y
N = 10; % Third dimension/depth of Y
Y = zeros(C,T,N);
%%
clc
n_movies = 2;
lastM = 0;
t = 0;
for p = 1:N
    m = ceil(p/(N/n_movies));
    if lastM == m
        % increment in time axis
        t = t + T;
        Y(:,:,p) = p1_data(m, 1:C, (t + 1):(t + T));        
    else
        t = 0;
        Y(:,:,p) = p1_data(m, 1:C, 1:T);
    end
    lastM = m;
end
%%
nice_plot(Y, 'Y3D', 0)
%%
K = 20; % Third dimension/depth of Y
D = zeros(C, T, K);
S = zeros(N, K);
gm_d = zeros(C, 1);
gm_s = 0;
gm_n = 0;

%% FIRST SAMPLES
clc
% Initialize gm_d
alpha_d = 3;
beta_d = 0.3;
for i = 1:C
    gm_d(i) = gamrnd(alpha_d, 1/beta_d);
end

% Initialize gm_s
alpha_s = 3;
beta_s = 0.3;
gm_s = gamrnd(alpha_s, 1/beta_s);

% Initialize gm_n
alpha_n = 3;
beta_n = 0.3;
gm_n = gamrnd(alpha_n, 1/beta_n);

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
        S(p, k) = normrnd(0, 1/sqrt(gm_s));
    end
end
%%
figure(1)
clf
imagesc(D(:,:,1));
%% INITIATE GIBBS
iters = 150;
burn = 75;
meanD = zeros(C,T,K);
meanS = zeros(N,K);
for it = 1:iters
    disp(['Iteration:' num2str(it)])
    if it > burn
        meanD = meanD + D;
        meanS = meanS + S;
    end
    [D, S, gm_d, gm_s, gm_n] = GibbsSampleNextTensor(Y, D, S, gm_d, gm_s, gm_n, constants, true);
end
%%
meanD = meanD/(iters - burn);
meanS = meanS/(iters - burn);
Y_approx = zeros(C, T, N);
for j = 1:T
    Y_approx(:,j,:) = squeeze(meanD(:,j,:))*meanS';
end
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
%%
figure(1)
clf
for k = 1:K
    imagesc(D(:,:,k)), colorbar
end
%%
figure(1)
clf
imagesc(S)










