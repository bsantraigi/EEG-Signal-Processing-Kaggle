%%
close all
clear
N = 10000;
X = normrnd(.5, 1, N, 1);

%%
alpha = 1e-3;
beta = 1e-3;

iter = 100000;
m = zeros(iter, 1);
gm = zeros(iter, 1);
m(1) = normrnd(0, 1);
gm(1) = gamrnd(alpha, 1/beta);

sum_x = sum(X);

for i = 2:iter
    gm_o = gm(i - 1);
    m(i) = normrnd(sum_x/(N + 1/gm_o), 1/sqrt(gm_o*N + 1));
    gm(i) = gamrnd(alpha + N/2, 1/(norm(X-m(i))^2/2 + beta));
end
sig = sqrt(1./gm);
%%
figure(1)
clf
histogram(m(2000:20:iter), 30)
title('Mu histogram');
figure(2)
clf
histogram(sig(2000:20:iter), 30)
title('var histogram');
