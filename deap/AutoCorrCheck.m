clearvars
%%
load('s05.mat')
%%
s = squeeze(data(1, 2, :));
hlag = 100;
xmat = zeros(2000, hlag); % Upto a delay of 100
for i = 1:hlag
    xmat(:, i) = s(i:(i + 1999));
end
%%
v = var(xmat(:, hlag));
c = cov(xmat);
%%
t = 0.8;
ps = t.^(0:(hlag-1));
figure(1)
clf
hold on
stem(flipud(c(:, hlag))/v)
plot(ps);
hold off
%%
load('s01.mat')
% Check signals from person 2 
% Validate the channel correlation