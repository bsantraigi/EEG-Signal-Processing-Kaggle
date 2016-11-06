function [ Y_new ] = repatch( muB, r, p, K, T )
%REPATCH2DATA Summary of this function goes here
% We had N = T*(r-p+1)^2 patches for all images combined
% Repatch will double the implicit patch size, p to 2*p and
% will generate the corresponding data matrix
% New size of data matrix would be N_new = T*(r-2*p+1)^2

% In muB each patch i.e. each column is represented by a 
% K dimensional vector
L = r - p + 1;
L_new = r - 2*p + 1;
Y_new = zeros(4*K, T*L_new^2);

for t = 1:T
    for u = 0:(L_new - 1)
        for v = 1:L_new
            j = v + u*L; % Previous patch index
            j_new = v + u*L_new; % New patch index
            newPatch = [muB(:, j), muB(:, j + p);...
                muB(:, j + p*L), muB(:, j + p + p*L)];
            Y_new(:, (t - 1)*L_new^2 + j_new) = newPatch(:);
        end
    end
end

end