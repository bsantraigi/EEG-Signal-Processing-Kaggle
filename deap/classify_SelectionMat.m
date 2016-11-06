clearvars
%%
load('Estimate_SMat_MultiUser[1234].mat');
%%
load('s01.mat')

[~, mGroups, uGroups] = MultiUserData([1], 32, 10, 40);
clear data
%%
valence = labels(:, 1);
arousal = labels(:, 2);
happy = find(valence > 6 & arousal > 5);
sad = find(valence < 6 & arousal < 5);
%%
figure(1)
clf
hold on
scatter(valence(sad), arousal(sad), 'r');
scatter(valence(happy), arousal(happy), 'g');
hold off
xlabel('Valence')
ylabel('Arousal')
%%
N = size(meanS, 1);
targets = zeros(N, 1);
for r = 1:N
    m = mGroups(r);
    if valence(m) > 6
        if arousal(m) > 6
            targets(r) = 1; % happy
        else
            targets(r) = 2; % relax
        end
    else
        if arousal(m) > 6
            targets(r) = 3; % angry
        else
            targets(r) = 4; % sad
        end
    end
end

t = [meanS, targets];