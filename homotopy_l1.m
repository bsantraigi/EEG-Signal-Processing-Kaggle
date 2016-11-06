%% 
clc
clear
close all
xs = normrnd(0, 1, 1000, 1);
xs= sort(xs);
system_order = 6;
a = [0;0;0.4;0;0;1];

X = ones(size(xs, 1),1);
for i = 1:system_order-1
    X = [X xs.^i];
end

y0 = X*a;
y = awgn(y0, 0.1, 'measured');
% y = y0;
mse2 = @(x, x2) mean((x-x2).^2);
plot(xs, y);


%%
clc
format long g
% j = 0
new_order = 8;

X_dat = ones(size(xs, 1), new_order);
for i = 1:new_order-1
    X_dat(:,i) = xs.^i;
end

%
a_j = zeros(new_order, 1);
lambda_j = max(abs(X_dat'*y));
fprintf('lambda_0:%.2f\n', lambda_j)
% Active elements contains same value as its index o/w zero
S_j = zeros(1, new_order); 

% FIND OUT THE FIRST ACTIVE ELEMENT
c_j = -X_dat'*y;
[m, l_1] = max(abs(c_j));

S_j(l_1) = l_1; % Add the new active element

% Find direction
% Find next breakpoint
% Update S_j, lambda_j
% Update a_j

% disp('------------------------------------')
toBreak = false;
iters = 4;
for j = 1:iters
%     S_j
%     c_j
    active_loci = (1:new_order == S_j);
    active_loci_c = (1:new_order ~= S_j);

    f_1 = X_dat(:, active_loci)'*X_dat(:, active_loci);
    d_j = zeros(new_order, 1);
    d_j(active_loci,1) = -f_1\sign(c_j(active_loci, 1));
    gm = zeros(new_order, 1);
    
    f_2 = X_dat'*X_dat*d_j;
    
    for l = 1:new_order
        if S_j(l) == l && j > 1
            % In active set
            v = -a_j(l, 1)/d_j(l, 1);
            if v <= 0
                v = Inf;
            end
        else
            % Not in active set
            
            v1 = ((lambda_j - c_j(l, 1))/(1 + f_2(l, 1)));
            v2 = ((lambda_j + c_j(l, 1))/(1 - f_2(l, 1)));
            if v1 < 0
                v1 = Inf;
            end
            if v2 < 0
                v2 = Inf;
            end
                
            v = min(v1, v2);
        end
        gm(l,1) = v;
    end
    [gm_j, new_i] = min(gm);
%     gm
%     fprintf('\ngm_j:%.2f, new_i:%.2f\n', gm_j, new_i);

    % Update S_j
    if S_j(new_i) > 0
%         disp({'Removing', new_i})
        S_j(new_i) = 0; % Remove element
    else
%         disp({'Adding', new_i})
        S_j(new_i) = new_i; % Add element
    end
        
    if gm_j > lambda_j
        disp({'BURN'})
        gm_j = lambda_j;
        toBreak = true;
    end
    % Update a_j
    a_j = a_j + gm_j*d_j;

    % Update c_j
    c_j = X_dat'*(X_dat*a_j - y);
    % Update lambda_j now (always at last), uses new a_j
    lambda_j_old = lambda_j;
    lambda_j = max(abs(c_j));

    fprintf('lambda_j:%.2f vs calculated:%.2f\n----------------------\n'...
        , lambda_j, lambda_j_old - gm_j)
%     if floor(lambda_j) ~= floor(lambda_j_old - gm_j)
%         disp('WHOA')
% %         break
%     end
    if toBreak
        break
    end
end

% c_j

figure(1)
clf
a_j
plot(xs, y, 'Color', [0.9,0.9,0.9]);
hold on;
plot(xs, y0,'blue');
y_approx = X_dat*a_j;
plot(xs, y_approx, 'red');

hold off;
legend('Noisy Samples', 'Ideal', 'Predicted');




