function alpha = Samplealpha(S,e0,f0,Z,alpha)
%Sample alpha
%Version 1: 09/12/2009
%Version 2: 10/21/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
ei = e0+0.5*numel(Z);
%fi = f0+0.5*sum(sum(S.^2));
fi = f0 + 0.5*sum(sum(S.^2)) + 0.5*(numel(Z)-nnz(Z))*(1/alpha);
alpha = gamrnd(ei,1./fi);
end