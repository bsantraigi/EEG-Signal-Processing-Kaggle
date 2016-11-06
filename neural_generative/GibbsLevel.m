function [ D, S, B, PI, post_PI, bias, Gamma ] = GibbsLevel( Y, D, S, B, PI, post_PI, bias, Gamma, Alpha, Beta, c )
%GIBBSLEVEL Summary of this function goes here
%   Detailed explanation goes here

D = sampleD(Y, D, S, B, bias, Gamma, c);
S = sampleS(Y, D, S, B, bias, Gamma, c);
[B, post_PI] = sampleB(Y, D, S, B, PI, post_PI, bias, Gamma, c);
PI = samplePI(B, PI, Alpha, Beta, c);
Gamma = sampleGammas(Y, D, S, B, bias, Gamma, Alpha, Beta, c);
bias = sampleBias(Y, D, S, B, Gamma, c);

end

