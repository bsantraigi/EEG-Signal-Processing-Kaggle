%%
close all
clear
N = 1000;
X = normrnd(.5, 1, N, 1);

s = sum(X);
s2 = sum(X.*X);
%%
clc
u = 200;
lambda = 0.001;
alpha_0 = 100;
beta_0 = 2;
beta = beta_0;
alpha = alpha_0;
for i = 1:5
    fprintf('\nITERATION: %d\n', i);
    fprintf('[u] Mean: %0.2f\tVariance: %0.4f\n', u, 1/lambda)
    fprintf('[gamma] Mean: %0.2f\tVariance: %0.4f\n', alpha/beta, alpha/beta^2)
    
    u = alpha*s/(N*alpha + beta);
    lambda = N*alpha/beta + 1;
    
    alpha = alpha + N/2;
    beta = beta + 0.5*(s2 + N/lambda + u*u*N - 2*u*s);
%     fprintf('u: %0.2f,\tlambda: %0.2f,\talpha: %0.2f\tbeta: %0.2f\n', u, lambda, alpha, beta);
    
end

% fprintf('[u] Mean: %0.2f\tVariance: %0.4f\n', u, 1/lambda)
% fprintf('[gamma] Mean: %0.2f\tVariance: %0.4f\n', alpha/beta, alpha/beta^2)
