%% 
clear
close all
N = 1000;
xs = normrnd(0, 5, N, 1);
xs = sort(xs);
system_order = 3;
a = [7;0.2;15]; % [a_0, a_1, a_2]
X = ones(size(xs, 1),1);
for i = 1:system_order-1
    X = [X xs.^i];
end
y0 = X*a;
% y = awgn(y0, 0.005, 'measured');
y = y0 + normrnd(0, 200, 1000, 1);
mse2 = @(x, x2) mean((x-x2).^2);

%%
clf
format long g
new_order = system_order; % Altered
% new_order = 8;
alpha_0 = 1e-1;
% beta_0 = 0.4; % Altered
beta_0 = 1e-1;
alpha_n = 1e-1;
beta_n = 1e-1;
gm_0 = gamrnd(alpha_0, 1/beta_0)
a_new = mvnrnd(zeros(new_order, 1), 1/gm_0*eye(new_order))
gm_n = gamrnd(alpha_n, 1/beta_n)

gsamples = 5000;

X_dat = X;

xtx = X_dat'*X_dat;

gm_0s = zeros(gsamples, 1);
as = zeros(new_order, gsamples);
gm_ns = zeros(gsamples, 1);

gm_0s(1) = gm_0;
gm_ns(1) = gm_n;
as(:, 1) = a_new;
% disp(as(:, 1)')
figure(1)
for it = 2:gsamples
%     disp(it)
    M = xtx + gm_0s(it - 1)/gm_ns(it - 1)*eye(new_order);
    u = M\X_dat'*y;
    sig = (1/gm_ns(it - 1))*inv(M);
    
    as(:, it) = mvnrnd(u, sig);
%     disp(as(:, it)')
    gm_0s(it) = gamrnd(alpha_0 + new_order/2, beta_0 + as(:, it)'*as(:, it)/2);
    
    gm_ns(it) = gamrnd(alpha_n + gsamples/2, beta_n + norm(y - X_dat*as(:, it))^2/2);
    
%     if it > 500
%         y_approx = X_dat*mean(as(:, 500:it), 2);
%         clf
%         hold on
%         plot(xs, y0);
%         plot(xs, y_approx);
%         hold off
%         drawnow
%     end
end
a_approx = mean(as(:,500:gsamples), 2)
clf
figure(1)
y_approx = X_dat*a_approx;
plot(xs, y);
hold on
plot(xs, y0);
plot(xs, y_approx);
legend('Noisy', 'Exact', 'Predicted Polynomial')

%%
clf
alpha = 1;
beta = 0.4;
figure(1)
z = gamrnd(alpha, 1/beta, 10000, 1);
mean(z)
var(z)
histogram(z)
