%%
user = 2;
f1 = load(sprintf('data/train_%d/%d_12_0.mat', user, user));
f1 = f1.dataStruct;
f2 = load(sprintf('data/train_%d/%d_10_1.mat', user, user));
f2 = f2.dataStruct;

%%
figure(2);
subplot(2,2,1)
plot(f1.data(:,1));
subplot(2,2,2)
f = abs(fft(f1.data(:,1)));
plot(f)
%%
a = 1;
b = 30;

d = 10240;
u = zeros(d, 1);
C = zeros(d, d);
b2 = b^2;
tic
for i = 1:d
    for j = 1:d
        C(i, j) = a*exp(-abs(i - j)/b2);
    end
end
toc
figure(1);
subplot(2, 2, 1);
imagesc(C);

y = mvnrnd(u, C);
subplot(2, 2, 2);
plot(y)

f = abs(fft(y));
subplot(2, 2, 3);
plot(f)

