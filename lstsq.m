%% 
clear
close all
xs = normrnd(0, 4, 500, 1);
xs = sort(xs);
a = [1.5;15;12];
X = [ones(size(xs, 1),1) xs xs.^2];
y0 = X*a;
y = awgn(y0, 1, 'measured');
mse2 = @(x, x2) mean((x-x2).^2);
%% LST-SQ
% load('lowrank_regsample.mat');
clc
figure(1)
clf
plot(xs, y0)
hold on 
plot(xs, y)
x_inv = (X'*X)\X';
disp('Least Square Result:');
a_approx = x_inv*y;
disp(a_approx')
y_approx = X*a_approx;
fprintf('MSE for LST Sq. %0.2f\n', mse2(y0, y_approx))
plot(xs, y_approx, 'r');
legend('Original', 'Noisy', 'Fitted')
title('Simple LST Sq')

% Ridge
figure(2)
clf
plot(xs, y0)
hold on 
plot(xs, y)
x_inv = (X'*X + 3800*eye(3))\X';
disp('Ridge Regression Result:');
a_approx = x_inv*y;
disp(a_approx')
y_approx = X*a_approx;
fprintf('MSE for RidgeR. %0.2f\n', mse2(y0, y_approx))
plot(xs, y_approx, 'r');
legend('Original', 'Noisy', 'Fitted')
title('With Ridge reg.')

%% Minimum RR

clc
figure(1)
clf
plot(xs, y0)
hold on 
plot(xs, y)
x_inv = (X'*X)\X';
disp('Least Square Result:');
a_approx = x_inv*y;
disp(a_approx')
y_approx = X*a_approx;
fprintf('MSE for LST Sq. %0.2f\n', mse2(y0, y_approx))
plot(xs, y_approx, 'r');
legend('Original', 'Noisy', 'Fitted')
title('Simple LST Sq')

% Ridge
figure(2)
clf
plot(xs, y0)
hold on 
plot(xs, y)
x_inv = (X'*X + 670*eye(3))\X';
disp('Ridge Regression Result:');
a_approx = x_inv*y;
disp(a_approx')
y_approx = X*a_approx;

last = 9999999;
for lam = 2000:50:4000
    x_inv = (X'*X + lam*eye(3))\X';
    a_approx = x_inv*y;
    y_approx = X*a_approx;
    le = mse2(y0, y_approx);
    if last > le
        last = le;
    else
        break
    end
    clf
    plot(xs, y0)
    hold on 
    plot(xs, y)
    plot(xs, y_approx, 'r');
    drawnow
end

fprintf('MSE for RidgeR. %0.2f (l = %0.1f)\n', le, lam)
legend('Original', 'Noisy', 'Fitted')
title('With Ridge reg.')