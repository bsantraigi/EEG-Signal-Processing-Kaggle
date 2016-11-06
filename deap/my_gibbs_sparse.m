clear
load('s01.mat')
%%
p1_data = data;
p1_labels = labels;

clear data labels
%%
n_movies = 2;
de = 32; % Number of channels
n_samples = size(p1_data, 3); % Num. of samples per movie
Y = zeros(de, n_movies*n_samples);
for i = 1:n_movies
    s = (i - 1)*n_samples + 1;
    f = i*n_samples;
    Y(:, s:f) = squeeze(p1_data(i, 1:32, :));
end
clear p1_data p1_labels
%%
frobNorm = @(M) sum(diag(M'*M));
K = 64;

n_iters = 4;
N = size(Y, 2);
D = zeros(de, K, n_iters);
S = zeros(K, N, n_iters); % Need to load-save simultaneously

gm_d = zeros(n_iters, 1);
alpha_d = 1;
beta_d = 1; % 0.01

gm_s = zeros(n_iters, 1);
alpha_s = 1;
beta_s = 1; % 0.003

gm_n = zeros(n_iters, 1);
alpha_n = 1; % 10
beta_n = 1;
% SAMPLE THE INITIAL VALUES

gm_d(1) = gamrnd(alpha_d, 1/beta_d);
gm_s(1) = gamrnd(alpha_s, 1/beta_s);
gm_n(1) = gamrnd(alpha_n, 1/beta_n);

D(:, :, 1) = normrnd(0, 1/sqrt(gm_d(1)), [de, K]);
S(:, :, 1) = normrnd(0, 1/sqrt(gm_d(1)), [K, N]);

for iter2 = 1:n_iters
    iter = 2; % Fake iter
    % Copy the results from this iteration into iter-1
    % Save all matrices and hyperparameters
    
    tic % START MEASURING TIME
    Q = S(:,:,iter - 1)*S(:,:,iter - 1)' + ...
        gm_d(iter - 1)/gm_n(iter - 1)*eye(K);
    D(:,:,iter) = ...
        matnormrnd(Y*S(:,:,iter - 1)'/Q, eye(de), inv(Q)/gm_n(iter - 1));
    Q = D(:,:,iter)'*D(:,:,iter) + gm_s(iter - 1)/gm_n(iter - 1)*eye(K);
    
    S(:,:,iter) = ...
        matnormrnd(Q\D(:,:,iter)'*Y, inv(Q)/gm_n(iter - 1), eye(N));
    gm_d(iter) = gamrnd(alpha_d + K*de/2, 1/(beta_d + ...
        frobNorm(D(:,:,iter))/2));
    gm_s(iter) = gamrnd(alpha_s + K*N/2, 1/(beta_s + ...
        frobNorm(S(:,:,iter))/2));
    gm_n(iter) = gamrnd(alpha_n + de*N/2, 1/(beta_n + ...
        frobNorm(Y - D(:,:,iter)*S(:,:,iter))/2));
    
    toc % REPORT EXECUTION TIME
    % Copy outputs back
    D(:,:,iter - 1) = D(:,:,iter - 1);
    S(:,:,iter - 1) = S(:,:,iter - 1);
    gm_n(iter - 1) = gm_n(iter);
    gm_d(iter - 1) = gm_d(iter);
    gm_s(iter - 1) = gm_s(iter);    
    
    % Saved the outputs to files
    filename = ['my_gibbs_matrices/iter_' num2str(iter2) '.mat'];
    IterResult = struct();
    IterResult.tempD = D(:,:,iter);
    IterResult.tempS = S(:,:,iter);
    IterResult.gmd = gm_d(iter);
    IterResult.gmn = gm_n(iter);
    IterResult.gms = gm_s(iter);
    save(filename, '-struct', 'IterResult');
    clear IterResult
    disp(['New file ' filename ' saved'])
end

%% READ SAVED FILES AND GET MEAN - WHAT TO DO ABOUT BURN-IN?
meanD = zeros(de, K);
meanS = zeros(K, N);
for iter2 = 2:n_iters
    filename = ['my_gibbs_matrices/iter_' num2str(iter2) '.mat'];
    load(filename)
    meanD = meanD + tempD;
    meanS = meanS + tempS;
    clear tempD tempS gmd gmn gms
end
meanD = meanD/n_iters;
meanS = meanS/n_iters;
%%

figure(1)
clf
Y_approx = meanD*meanS;
for c = 1:5
    plot(Y(:,c))
    hold on
    plot(Y_approx(:,c))
    hold off
    drawnow
    pause(1)
end
clear Y_approx

%%
figure(1)
clf
plot(meanD(:,1:3))