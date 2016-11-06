%%

% X = [1 2 0 0 0;
%     0 0 3 5 1;
%     4 5 0 1 2;
%     7 9 0 3 3;
%     4 5 2 1 0;
%     3 3 4 1 9;
%     4 0 4 0 4];
% x = - (cos(2*pi*2*t)-1).*exp(-5*t);

% X = X'; % 8 4D training vectors

D1 = rand(5, 5);
X = rand(5, 100);
%% Iterative Least Squares
D = D1;
for it = 1:1000
    W = sparseapprox(X, D, 'ORMP', 'tnz', 4);
    D = (X*W')/(W*W');
    D = dictnormalize(D);
end

W = sparseapprox(X, D, 'ORMP', 'tnz', 2);
disp(W)

%% K - SVD
D = D1;
for it = 1:10
    W = sparseapprox(X, D, 'ORMP', 'tnz',1);
    R = X - D*W;
    for k=1:2
        I = find(W(k,:));
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        % U is normalized
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end    
end

W = sparseapprox(X, D, 'ORMP', 'tnz', 1);
disp(W)