function [ X ] = RemoveComplex( X )
%REMOVECOMPLEX Summary of this function goes here
%   Columns are features

% last column is label - no fix required there
for f = 1:(size(X, 2))
    r = (X(:, f) == real(X(:,f)));
    m = mean(X(r, f));
    X(~r, f) = m;
end

