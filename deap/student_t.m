clear
clc
iters = 10000;
x_samples = zeros(iters, 1);
gm_samples = zeros(iters, 1);

% alpha_d = 0.625000;
% beta_d = 0.001250;

alpha_d = 200;
beta_d = 500000;

x_samples(1) = normrnd(0,1);
gm_samples(1) = gamrnd(alpha_d, 1/beta_d);
for j=2:iters
    gm_samples(j) = gamrnd(alpha_d + 0.5, 1/(beta_d + 0.5*x_samples(j-1)^2));
    x_samples(j) = normrnd(0, sqrt(1/gm_samples(j)));
    
%     gm_samples(j) = gamrnd(alpha_d, 1/beta_d);
%     x_samples(j) = normrnd(0, sqrt(1/gm_samples(j)));
end
figure(2)
clf
h = hist(x_samples(2000:size(x_samples, 1)), 100);
hist(x_samples(2000:size(x_samples, 1)), 100)
m = mean(x_samples(2000:size(x_samples, 1)))
v = var(x_samples(2000:size(x_samples, 1)))
hold on
z = normrnd(m, sqrt(v), 18000, 1);
% hist(z)
% axis([-5 5 0 max(h)*1.1])

%%
c = 1.5; % Decrease to increase variance
k = 0.5; % Increase to increase mean
alpha_d = c*k^2;
beta_d = 0.001*c*k;
fprintf('HI\n\nalpha_d = %f;\nbeta_d = %f;\n', alpha_d, beta_d)
z = gamrnd(alpha_d, 1/beta_d, 100000, 1);
fprintf('\nmean = %f;\nvar = %f;\n\n', mean(z), var(z))
fprintf('\nmean = %f;\nvar = %f;\n\n', alpha_d/beta_d, alpha_d/beta_d^2)
figure(2)
clf
hist(z, 100)
%%
c = 10; % Decrease to increase variance
k = 0.2; % Increase to increase mean
alpha_d = c*k^2;
beta_d = 0.0001*c*k;
fprintf('HI\n\nalpha_d = %f;\nbeta_d = %f;\n', alpha_d, beta_d)
z = gamrnd(alpha_d, 1/beta_d, 100000, 1);
fprintf('\nmean = %f;\nvar = %f;\n\n', mean(z), var(z))
fprintf('\nmean = %f;\nvar = %f;\n\n', alpha_d/beta_d, alpha_d/beta_d^2)
figure(2)
clf
hist(z, 100)
%%
c = 0.07; % Decrease to increase variance
k = 10; % Increase to increase mean
alpha_d = c*k^2;
beta_d = c*k;
fprintf('HI\n\nalpha_d = %0.10f;\nbeta_d = %0.10f;\n', alpha_d, beta_d)
z = gamrnd(alpha_d, 1/beta_d, 1000000, 1);
fprintf('\nmean = %f;\nvar = %f;\n\n', mean(z), var(z))
fprintf('\nmean = %f;\nvar = %f;\n\n', alpha_d/beta_d, alpha_d/beta_d^2)
figure(2)
clf
hist(z, 100)
%%
igamrnd = @(A, B) 1./gamrnd(A, B, 10000, 1);
zs0 = igamrnd(31, 1/16);

figure(1)
clf
% histogram(zs, 50, 'BinLimits', [0, 3]);
step = 0.2;
xlim = 2;
h0 = histc(zs0, (0:step:xlim));
plot(0:step:xlim, h0)

%%
zs0 = gamrnd(4, 1/2, 10000,1);

figure(1)
clf
% histogram(zs, 50, 'BinLimits', [0, 3]);
step = 0.2;
xlim = 20;

h0 = histc(zs0, (0:step:xlim));
plot(0:step:xlim, h0)
%%
aa = 31;
modee = 0.5;
stdd = 0.1;
ratiosq = (aa - 1)^2*(aa - 2)/(aa + 1)^2;
betta = (aa + 1)*modee;
fprintf('rtarget: %f, robtained: %f, beta: %f\n', (modee/stdd)^2, ratiosq, betta);

