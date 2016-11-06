%% Synthetic Dataset
clear
close all
t2 = 0:499;
t3 = 500:599;
X1 = [square(2*pi*t2/100), sin(2*pi*t2/100)]';
X2 = [sin(2*pi*t2/100), square(2*pi*t2/100)]';
X3 = [cos(2*pi*t2/100), sawtooth(2*pi*t2/100)]';
X4 = [sawtooth(2*pi*t2/100), cos(2*pi*t2/100)]';
figure()

subplot(2,2,1), plot(X1)
subplot(2,2,2), plot(X2)
subplot(2,2,3), plot(X3)
subplot(2,2,4), plot(X4)
title('Original Basis')
%%
sig = zeros(1000, 30);
sc = 5;

for i = 1:sc
    sig(:,i) = X1 + 0.25*(rand(1000,1) - 0.5);
end

for i = 1:sc
    sig(:,sc + i) = X2 + 0.25*(rand(1000,1) - 0.5);
end

for i = 1:sc
    sig(:,2*sc + i) = X3 + 0.6*(rand(1000,1) - 0.5);
end

for i = 1:sc
    sig(:,3*sc + i) = X4 + 0.2*(rand(1000,1) - 0.5);
end

for i = 1:sc
    sig(:,4*sc + i) = X4 + X1 + 0.77*(rand(1000,1) - 0.5);
end

for i = 1:sc
    sig(:,5*sc + i) = X2 + X3 + 0.3*(rand(1000,1) - 0.5);
end

figure()
subplot(3,2,1), plot(sig(:,1))
subplot(3,2,2), plot(sig(:,6))
subplot(3,2,3), plot(sig(:,12))
subplot(3,2,4), plot(sig(:,18))
subplot(3,2,5), plot(sig(:,22))
subplot(3,2,6), plot(sig(:,27))
title('Training Signals')
%% Real EEG Dataset
close all
clear
sig = load('deapSmallOrdered.mat');
sig = sig.X;

%%
%[D, Z, S]=BPFA_SDL(sig, maxIt, K, displayPeriod);
[D, Z, S] = BPFA_SDL(sig, 100, 8, 40);
figure()
imagesc(Z.*S)
%%
clf
figure(gcf)
imagesc(Z.*S)
%%
figure()
n = 3;
m = ceil(size(D,2)/n);

for j=1:size(D,2)
    subplot(m,n,j);
    plot(D(:,j))
end
%%
figure()
X_a = D*(Z.*S)';
for i=1:6
    subplot(3,2,i)
    plot(X_a(:,5*(i-1) + 3))
    hold on
    plot(sig(:,5*(i-1) + 3));
end