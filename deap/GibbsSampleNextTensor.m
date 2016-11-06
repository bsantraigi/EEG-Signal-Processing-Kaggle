function [D, S, gm_d, gm_s, gm_n] = GibbsSampleNextTensor...
    (Y, D, S, gm_d, gm_s, gm_n, constants, updateD )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [C, T, K] = size(D);
    [~, ~, N] = size(Y);
    
    alpha_d = constants.alpha_d;
    beta_d = constants.beta_d;
    alpha_s = constants.alpha_s;
    beta_s = constants.beta_s;
    alpha_n = constants.alpha_n;
    beta_n = constants.beta_n;
    
    if updateD
        %Sample D
        for j = 1:T
            for i = 1:C
                for k = 1:K
                    delY = squeeze(Y(i,j,:)) -...
                        S(:, [1:(k-1), (k+1):K])*squeeze(D(i, j, [1:(k-1), (k+1):K]));
                    if j == 1
                        numer = gm_d(i)*D(i, j+1, k) + ...
                            gm_n*sum(S(:,k).*(delY));
                    elseif j < T
                        numer = gm_d(i)*(D(i, j+1, k) + D(i, j-1, k)) + ...
                            gm_n*sum(S(:,k).*(delY));
                    else
                        numer = gm_d(i)*D(i, j-1, k) + ...
                            gm_n*sum(S(:,k).*(delY));
                    end
                    denom = 2*gm_d(i) + gm_n*sum(S(:, k).^2);
                    MU = numer/denom;
                    Sigma = 1/sqrt(denom);
                    D(i, j, k) = normrnd(MU, Sigma);
                end
            end
        end
    end
    
    % Sample S
    for k = 1:K
        Y_except_k = zeros(C, T, N);
        for j = 1:T
            Y_except_k(:,j,:) = ...
                squeeze(D(:,j,[1:(k - 1), (k + 1):K]))*S(:, [1:(k - 1), (k + 1):K])';
        end
        
        del_y = Y - Y_except_k;
        for p = 1:N        
            numer = gm_n*sum(sum(D(:,:,k).*del_y(:,:,p)));
            denom = norm(D(:,:,k), 'fro')^2 + gm_s;
            MU = numer/denom;
            Sigma = 1/sqrt(denom);
            S(p,k) = normrnd(MU, Sigma);
        end
    end
    
    % Sample gm_d
    del_d = D;
    del_d(:, 2:T, :) = del_d(:, 2:T, :) - del_d(:, 1:(T - 1), :);
    ALPHA = K*T/2 + alpha_d;
    
    for i = 1:C
        BETA = beta_d + 0.5*norm(squeeze(del_d(i, :, :)) ,'fro').^2;
        gm_d(i) = gamrnd(ALPHA, 1/BETA);
    end
    
    % Sample gm_s
    ALPHA = K*N/2 + alpha_s;
    BETA = beta_s + 0.5*norm(S, 'fro')^2;
    gm_s = gamrnd(ALPHA, 1/BETA);
    
    % Sample gm_n
    Y_approx = zeros(C, T, N);
    for j = 1:T
        Y_approx(:,j,:) = squeeze(D(:,j,:))*S';
    end
    ALPHA = C*N*T/2 + alpha_n;
    del_tensor = sum(sum(sum((Y - Y_approx).^2)));
    BETA = beta_n + 0.5 * del_tensor;
    gm_n = gamrnd(ALPHA, 1/BETA);

end

