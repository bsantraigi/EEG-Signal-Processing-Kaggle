%%
clear
X = load('s01.mat');
sig = X.data;

%% Plot raw data of nearby electrodes
% Found - Positive spatial correlations
close all
for i = 1:1000
    sig_partial = squeeze(sig(1,1:32,(1 + i):(100 + i)));
%     plot(sig_partial')
    imagesc(sig_partial')
    drawnow
end