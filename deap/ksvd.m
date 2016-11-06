function [ D, W ] = ksvd( X, dictSize, iter, tnz, reduceSize )
%KSVD Summary of this function goes here
%   Detailed explanation goes here

D = rand(size(X, 1), dictSize);
norm2 = sqrt(sum(D.*D));
D = D./repmat(norm2, size(X, 1), 1);

mag = @(v) sqrt(sum(v.^2, 1));

t = 1:iter;
errors = zeros(iter, 1);
if nargin < 4
    tnz = 5;
end
if nargin < 5
    reduceSize = false;
end
% figure(gcf);
%% K - SVD original
for it = 1:iter
    W = sparseapprox(X, D, 'OMP', 'tnz',tnz);
    X_a = D*W;
    R = X - X_a;
    if(reduceSize)
        useless = find(sum(W~=0, 2) == 0);
        D(:, useless) = [];
        W(useless, :) = [];
        dictSize = size(D, 2);
        disp(['Reduced dictionary size ' , num2str(dictSize)])
    end
%     clf;
%     plot(X(3000:4000,181),'g')
%     hold on
%     plot(X_a(3000:4000, 181),'r')
%     drawnow;
    
    for k=1:dictSize
        I = find(W(k,:));
        if size(I, 2) == 0
            continue
        end
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        % U is normalized
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end
%     er = mean(mag(R));
%     disp(['Error ', num2str(it), ':', num2str(er)])
%     errors(it) = mean(mag(R));
    snrs = snr_2(X, D*W);
    fprintf('Range of SNR %d: (%f, %f)\n', it, min(snrs), max(snrs));
    
%     plot(t, errors);
%     drawnow;
    
end
W = sparseapprox(X, D, 'OMP', 'tnz',tnz);

end

