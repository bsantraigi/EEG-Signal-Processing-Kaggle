clearvars
load('s02.mat');
%%
figure(1)
    clf
for frame = 1:4
    Y = squeeze(data(frame, 1:32, 1:512));
    subplot(2, 2, frame)
    imagesc(Y), colorbar
end
%%
fz = @(x) 1/det(GetCovMat(x, 0.2, T));
figure(1)
clf
hold on
for e = 0.99:0.01:1.1
    plot(e, fz(e), 'x')
end
hold off
%% Variable Defs.
% Inverse gamma dist.
igamrnd = @(A, B) 1/gamrnd(A, B);

T = 512;
C = 32;
K = 50;
Y = squeeze(data(1, 1:C, 1:T))';

Alpha = struct();
Beta = struct();

D = zeros(T, K);
AllSigmas = zeros(T, T, K); % K covariance matrices
ijMat = (repmat(1:T, T, 1) - repmat((1:T)', 1, T)).^2;
GetCovMat = @(ak, Bk, T) ak.*exp(-ijMat./Bk);

aMat = zeros(K, 1);
Alpha.a = 10;
Beta.a = 5;

BMat = zeros(K, 1);
Alpha.b = 4;
Beta.b = 2;

% Initial Values of a
% Initial Values of B
% Initial Values of Sigma
for k = 1:K 
    aMat(k) = igamrnd(Alpha.a, 1/Beta.a);
    BMat(k) = gamrnd(Alpha.b, 1/Beta.b);
    AllSigmas(:,:,k) = GetCovMat(aMat(k), BMat(k), T);
%     AllSigmas(:,:,k) = GetCovMat(10, 0.2, T);
    D(:, k) = mvnrnd(zeros(T, 1), AllSigmas(:,:,k));
    mvnvalue = mvnpdf(D(:, k), zeros(T, 1), AllSigmas(:,:,k));
    fprintf('[%d] %f, %f, %f, %f \n', k, aMat(k), BMat(k), AllSigmas(1,1,k), mvnvalue);
end

Alpha.n = 10;
Beta.n = 4;
Gamma.n = gamrnd(Alpha.n, 1/Beta.n);

S = zeros(K, C);
Alpha.s = 1.5;
Beta.s = 0.005;
Gamma.s = zeros(K, C);

for k = 1:K
    for c = 1:C
        Gamma.s(k, c) = gamrnd(Alpha.s, 1/Beta.s);
        S(k, c) = normrnd(0, 1/sqrt(Gamma.s(k, c)));
    end
end

%% Start Sampling
% [D, AllSigmas, aMat, BMat, S, Gamma, Alpha, Beta] ...
%     = GibbsSampleNext_GP_v1(Y, D, AllSigmas, aMat, BMat, S, Gamma, Alpha, Beta, true);
updateD = true;

% Sample D
I = eye(T);
if updateD
    lambda_d2 = Gamma.n*sum(S.^2, 2);
    for k = 1:K
        delY = Y - D(:, [1:(k - 1) (k+1):K])*S([1:(k - 1) (k+1):K], :);
        mu_dk2 = sum(repmat(S(k,:), T, 1).*delY(:,:), 2)/sum(S(k,:).^2);
        ldk2 = lambda_d2(k)*I;
%         det(AllSigmas(:,:,k))
        SIGD = inv(ldk2 + inv(AllSigmas(:,:,k))); % COV Mat
        MUD = SIGD*ldk2*mu_dk2;
        D(:,k) = mvnrnd(MUD, SIGD);
    end
end
%% Sample a

for k = 1:K
    % P.S. - The a_k factor has been cacelled by multiplication
    beta_temp = Beta.a + 0.5*aMat(k)*(D(:, k)'/AllSigmas(:,:,k))*D(:,k);
    aMat(k) = igamrnd(Alpha.a + T/2, 1/beta_temp);
end
%% Sample B
bpdf = @(next) mvnpdf(D(:,k), zeros(T, 1), GetCovMat(aMat(k), next, T))*gampdf(next, Alpha.b, 1/Beta.b);
bprop = @(last, next) gampdf(next, Alpha.b, 1/Beta.b);
bproprnd = @(last) gamrnd(Alpha.b, 1/Beta.b);
for k = 1:K
    b_samples = mhsample(BMat(k), 500, 'pdf', bpdf, 'proppdf', bprop, 'proprnd', bproprnd);
    BMat(k) = mean(b_samples(200:500));
end

%% Sample S
for k = 1:K
    delY = Y - D(:, [1:(k - 1) (k+1):K])*S([1:(k - 1) (k+1):K], :);
    H = sum(repmat(D(:,k), 1, C).*delY, 1);
    DsqSum = sum(D.^2, 1);
    for c = 1:C
        VARS = 1/(Gamma.n*DsqSum(k) + Gamma.s(k,c)); % Std. Dev.
        MUS = Gamma.n*H(c)*VARS;
        S(k,c) = normrnd(MUS, sqrt(VARS));
    end
end

% Sample Gamma.s
for k = 1:K
    for c = 1:C
        Gamma.s(k, c) = gamrnd(Alpha.s + 0.5, 0.5*S(k,c)^2 + Beta.s);
    end
end

% Sample Gamma.n
delY = sumsqr(Y-D*S);
Gamma.n = gamrnd(T*C/2 + Alpha.n, Beta.n + 0.5*delY);
%%
figure(1)
clf
imagesc(S), colorbar
%%

Y_approx = D*S;
figure(2)
subplot(1,2,1)
imagesc(Y)
subplot(1,2,2)
imagesc(Y_approx)
corr(Y(:), Y_approx(:))
%%
figure(2)
subplot(1,2,1)
imagesc(Y)
subplot(1,2,2)
imagesc(D)
%%
figure(1)
clf
hold on
plot(Y(:,2))
plot(D(:,10))
hold off
legend('Y', 'D');