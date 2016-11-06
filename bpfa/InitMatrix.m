function [D,S,Z,phi,alpha,Pi] = InitMatrix(X,P,N,K,InitOption,UpdateOption, Yflag,IsSubMean,MuX1)
%Initialization
%Version 1: 09/12/2009
%Version 2: 10/21/2009
%Version 3: 10/26/2009
%Version 4: 10/28/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
if nargin<7
    Yflag=[];
    IsSubMean = [];
    MuX1 = [];
else
    if ~IsSubMean
        X = X.*Yflag + repmat(MuX1,1,N).*(~Yflag);
    end
end
phi = 1/((25/255)^2);
alpha = 1;
if strcmp(InitOption,'SVD')==1 || strcmp(InitOption,'SVD1')==1
    [U_1,S_1,V_1] = svd(full(X),'econ');
    if P<=K
        D = zeros(P,K);
        D(:,1:P) = U_1*S_1;
        S = zeros(N,K);
        S(:,1:P) = V_1;
    else
        D =  U_1*S_1;
        D = D(1:P,1:K);
        S = V_1;
        S = S(1:N,1:K);
    end
    if strcmp(InitOption,'SVD')==1
        Z = true(N,K);
        Pi = 0.5*ones(K,1);
    else
        Z = sparse(false(N,K));
        Pi = 0.01*ones(K,1);
    end
elseif strcmp(InitOption,'RAND')==1
    D=randn(P,K)/sqrt(P);
    S=randn(N,K);
    Z = sparse(false(N,K));
    Pi = 0.01*ones(K,1);
elseif  strcmp(InitOption,'Kmeans')==1
    S = ones(N,K);
    Z = sparse(false(N,K));
    if N<K;
        Xtemp=zeros(P,K);
        Xtemp(:,1:N)=X;
    else
        Xtemp=X;
    end
    [idx,ctrs] = kmeans(Xtemp',K,'emptyaction','singleton');
    D  = ctrs';
    Pi = 0.01*ones(K,1);
    if nargin<7
        [X_k, D, Z, S] = SampleDZS(X - D*(Z.*S)', D, Z, S, Pi, alpha, phi, K, false, true, true,UpdateOption);
    else
        [X_k, D, Z, S] = SampleDZS_MissingData( (X - D*(Z.*S)').*Yflag, D, Z, S, Pi, alpha, phi, P, K, N, Yflag, false, true, true,UpdateOption);
    end
end
end





