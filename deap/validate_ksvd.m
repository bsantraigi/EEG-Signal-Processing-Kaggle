%% Validate KSVD
clear

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
S = zeros(1000, 20);
sc = 5;

for i = 1:sc
    S(:,i) = X1 + 0.4*(rand(1000,1) - 0.5);
end

for i = 1:sc
    S(:,sc + i) = X2 + 0.4*(rand(1000,1) - 0.5);
end

for i = 1:sc
    S(:,2*sc + i) = X3 + 0.4*(rand(1000,1) - 0.5);
end

for i = 1:sc
    S(:,3*sc + i) = X4 + 0.4*(rand(1000,1) - 0.5);
end

figure()
subplot(2,2,1), plot(S(:,1))
subplot(2,2,2), plot(S(:,6))
subplot(2,2,3), plot(S(:,12))
subplot(2,2,4), plot(S(:,18))
title('Training Signals')
%%
[u,s,v] = svd(S);
s = max(s);
% fprintf('Singular value range of X: %f, %f\n', min(s), max(s));
figure()
plot(s)
title('Singular Values')
%%
[D, W] = ksvd(S, 6, 50);
%%
figure()
for j=1:size(D,2)
    subplot(3,2,j);
    plot(D(:,j))
end
title('Dictionary Atoms')
%%
S_a = D*W;
figure()
subplot(2,2,1), plot(S(:,1))
hold on, plot(S_a(:,1))
subplot(2,2,2), plot(S(:,6))
hold on, plot(S_a(:,6))
subplot(2,2,3), plot(S(:,12))
hold on, plot(S_a(:,12))
subplot(2,2,4), plot(S(:,18))
hold on, plot(S_a(:,18))
title('Reconstructed Signals')