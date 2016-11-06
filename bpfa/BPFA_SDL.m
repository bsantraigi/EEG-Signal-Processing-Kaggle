function [ D, Z, S ] = BPFA_SDL( X_k, maxIt, K, displayPeriod )
%BPFA_SDL Learn a dictionary for sparse representation of the data
%   Detailed explanation goes here
% Manually set variables

[P,N] = size(X_k);
if nargin < 3
    K = N;
end
if nargin < 4
    displayPeriod = maxIt/5;
end
InitOption = 'SVD';
UpdateOption = 'DkZkSk';

%Sparsity Priors        
if strcmp(InitOption,'SVD')==1
    a0 = 1;
    b0 = N/8;
else
    a0=1;
    b0=1;
end
% Set Hyperparameters
c0=1e-6;
d0=1e-6;
e0=1e-6;
f0=1e-6;

NoiseVar = [];
ReduceDictSize = true;

%Initializations for new added patches
[D,S,Z,phi,alpha,Pi] = InitMatrix(X_k,P,N,K,InitOption,UpdateOption);

X_k = X_k - D*(Z.*S)';
for iter=1:maxIt
    tic
    [X_k, D, Z, S] = SampleDZS(X_k, D, Z, S, Pi, alpha, phi, P, K, N, true, true, true, UpdateOption);
    Pi = SamplePi(Z,N,a0,b0);
    alpha = Samplealpha(S,e0,f0,Z,alpha);
    phi = Samplephi(X_k,c0,d0);
    ittime=toc;
    if(mod(iter, displayPeriod) == 0)
        disp(['Chekpoint:', num2str(iter), '; ', ...
            'Dictionary Size:', num2str(K)])
    end

    NoiseVar(end+1) = sqrt(1/phi)*255;
    if ReduceDictSize
        sumZ = sum(Z,1)';
        if min(sumZ)==0
            Pidex = sumZ==0;
            D(:,Pidex)=[];
            K=size(D,2);
            S(:,Pidex)=[];
            Z(:,Pidex)=[];
            Pi(Pidex)=[];
        end
    end    
end

end

