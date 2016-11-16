function [ L, logLike0, logLike1 ] = CalcLogLR( fs, histMat_0, histEdges_0, histMat_1, histEdges_1 )
%CALCLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
L = 0;
histEdges_0 = histEdges_0(:, 1:end-1);
histEdges_1 = histEdges_1(:, 1:end-1);
m1 = min(histMat_1, 2);
m0 = min(histMat_0, 2);
logLike1 = 0;
logLike0 = 0;
for f = 1:256
    k0 = interp1(histEdges_0(f, :), histMat_0(f, :), fs(f));
    k1 = interp1(histEdges_1(f, :), histMat_1(f, :), fs(f));
    if k1 == 0
        k1 = m1(f);
    end
    if k0 == 0
        k0 = m0(f);
    end
    if isnan(k1)
        if isnan(k0)
            k1 = 1e-2;
            k0 = 1e-2;
        else
            k1 = k0/2;
        end
    end
    if isnan(k0)
        if isnan(k1)
            k0 = 1e-2;
            k1 = 1e-2;
        else
            k0 = k1/2;
        end
    end
    logLike0 = log(k0) + logLike0;
    logLike1 = log(k1) + logLike1;
%     fprintf('Inner L: %f\n', L);
    L = L + log(k1/k0);
end
end

