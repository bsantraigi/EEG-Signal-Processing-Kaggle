%% Actual Code
clear all
close all
imgPath = 'caltech101/101_ObjectCategories/'
typeofimage = 'Faces_easy/'
fl = dir([imgPath typeofimage]);
%%
figure(1)
Y = [];
reduceTo = 32;
patchsize = 4;
totalImages = 2;
gap = 100;
for imindex = 3:gap:(3 + gap*totalImages - 1)    
    imgTemp = imread([imgPath typeofimage fl(imindex).name]);
    if size(imgTemp, 3) > 1
        imgTemp = rgb2gray(imgTemp);
    end
    imgTemp = imresize(imgTemp, reduceTo/min(size(imgTemp)));
    if size(imgTemp, 1) > reduceTo
        r = size(imgTemp, 1);
        drop = floor((r - reduceTo)/2);
        s = max(1, 1+drop - 1);
        imgTemp = imgTemp(s:(s+reduceTo-1),:);
    elseif size(imgTemp, 2) > reduceTo
        r = size(imgTemp, 2);
        drop = floor((r - reduceTo)/2);
        s = max(1, 1+drop - 1);
        imgTemp = imgTemp(:,s:(s+reduceTo-1));
    end
    Y = [Y im2patch(imgTemp, patchsize)];
    size(imgTemp)
    imshow(imgTemp)
    drawnow
    pause(0.3)
end
Y = Y./255;
clearvars r s drop
%%
figure(1)
clf

recon = patch2im(Y(:,1:size(Y,2)/totalImages), patchsize);
% recon = reshape(D(:, 23), reduceTo, reduceTo);
imshow(recon)
%% Initialize
K1 = 40;
K2 = 60;

Alpha1 = {};
Beta1 = {};

% Params for gamma distro
Alpha1.d = 1e-1;
Beta1.d = 1e-1;
Alpha1.s = 1e-1;
Beta1.s = 1e-1;
Alpha1.bias = 1e-1;
Beta1.bias = 1e-1;
Alpha1.n = 1e-3;
Beta1.n = 1e-3;

% Params for beta distro : Near to zero, sparse
Alpha1.pi = 1;
Beta1.pi = 1000;

Alpha2 = Alpha1;
Beta2 = Beta1;

Alpha2.pi = 1;
Beta2.pi = 2000;

pi2y = @(PI_Mat) log(PI_Mat./(1-PI_Mat));

[ D, S, B, PI, post_PI, bias, Gamma, c ] = InitAll( Y, K1, Alpha1, Beta1 );

Y2 = sigmoid_Inv(repatch(post_PI, reduceTo, patchsize, K1, totalImages));
[ D2, S2, B2, PI2, post_PI2, bias2, Gamma2, c2 ] = InitAll( Y2, K2, Alpha2, Beta2 );

%% Gibbs - Only Layer 1
figure(2)
clf

tune_length = 100;
round_1 = tune_length;
round_2 = round_1 + tune_length;

for gr = 1:2000
    % Test here only
    Y_approx = D*(S.*B) + repmat(bias, 1, c.N);
    l = (reduceTo - patchsize + 1)^2;
    for r = 0:l:1
        
        subplot(3, 2, 1)
        recon = patch2im(Y_approx(:,(r+1):(r+l)), patchsize);
        imshow(recon)
        title('Recon')
        
        subplot(3, 2, 2)
        imshow(patch2im(Y(:,(r+1):(r+l)), patchsize))
        title('Actual')
        
        subplot(3, 2, 3)
        imagesc(B), colorbar
        title('B Matrix')
        
        subplot(3, 2, 4)
        imagesc(Y2), colorbar
        title('Y2 Matrix')
        
        subplot(3, 2, 5)
        imagesc(B2), colorbar
        title('B2 Matrix')
        
        subplot(3, 2, 6)
        imagesc(post_PI2), colorbar
        title('PI2 Matrix')
        drawnow
    end
    tic

    if gr < round_1
        % LEarn layer 1
        [ D, S, B, PI, post_PI, bias, Gamma] = GibbsLevel( Y, D, S, B, PI, post_PI, bias, Gamma, Alpha1, Beta1, c );
        Y2 = sigmoid_Inv(repatch(post_PI, reduceTo, patchsize, K1, totalImages));
    elseif gr < round_2
        %Learn Layer 2
        [ D2, S2, B2, PI2, post_PI2, bias2, Gamma2] = GibbsLevel( Y2, D2, S2, B2, PI2, post_PI2, bias2, Gamma2, Alpha2, Beta2, c2 );
    else
        % Learn Both
        [ D, S, B, PI, post_PI, bias, Gamma] = GibbsLevel( Y, D, S, B, PI, post_PI, bias, Gamma, Alpha1, Beta1, c );
        Y2 = sigmoid_Inv(repatch(post_PI, reduceTo, patchsize, K1, totalImages));
        [ D2, S2, B2, PI2, post_PI2, bias2, Gamma2] = GibbsLevel( Y2, D2, S2, B2, PI2, post_PI2, bias2, Gamma2, Alpha2, Beta2, c2 );
    end
    
    % Checkpoint for B - Layer 1
    if mod(gr, 2) == 0
        if sum(sum(B == 0)) == c.N*K1
            display('Resetting B1')
            [ ~, ~, B, ~, ~, ~, ~, ~ ] = InitAll( Y, K1, Alpha1, Beta1 );
        end
    end
    
    % Layer 2 Checkpoint
    if mod(gr, 5) == 0
        if sum(sum(B2 == 0)) == c.N*K2
            display('Resetting B2')
            [ ~, ~, B2, ~, ~, ~, ~, ~ ] = InitAll( Y2, K2, Alpha2, Beta2 );
        end
    end
    
    
    toc
    fprintf('Iteration: %d \n', gr)
    fprintf('Noise Var: L1 -> %f, L2 -> %f\n', 1/Gamma.n, 1/Gamma2.n)
end
fprintf('Gibbs Complete...\n')

%% Plot reconstructed Image - Using Layer 1
figure(2)
clf
Y_approx = D*(S.*B) + repmat(bias, 1, c.N);
l = (reduceTo - patchsize + 1)^2;
for r = 0:l:2500
    subplot(1, 2, 1)
    recon = patch2im(Y_approx(:,(r+1):(r+l)), patchsize);
    recon(recon<=0) = 0;
    recon(recon>=1) = 1;
    imshow(recon)
    title('Recon')
    subplot(1, 2, 2)
    imshow(patch2im(Y(:,(r+1):(r+l)), patchsize))
    title('Actual')
    drawnow
    pause(0.3)
end
%% Plot reconstructed Image - Layer 2
figure(2)
clf
Y2_approx = D2*(S2.*B2) + repmat(bias2, 1, c2.N);
muB_approx = un_repatch(1./(1+exp(-Y2_approx)), reduceTo, patchsize, K1, totalImages);
B_approx = binornd(ones(K1, c.N), muB_approx);
Y_approx = D*(S.*muB_approx) + repmat(bias, 1, c.N);
l = (reduceTo - patchsize + 1)^2;
for r = 0:l:2500
    subplot(1, 2, 1)
    recon = patch2im(Y_approx(:,(r+1):(r+l)), patchsize);
    recon(recon<=0) = 0;
    recon(recon>=1) = 1;
    imshow(recon)
    title('Recon')
    subplot(1, 2, 2)
    imshow(patch2im(Y(:,(r+1):(r+l)), patchsize))
    title('Actual')
    drawnow
    pause(0.3)
end
%% Plot the sorted B matrix
[~, bi] = sort(sum(B, 2));
figure(5), imagesc(B(bi, :))
%% Plot the sorted B2 matrix
[~, bi] = sort(sum(B2, 2));
figure(5), imagesc(B2(bi, :))
%% Plot all features
normalize = @(mat) (mat - min(min(mat)))/(max(max(mat)) - min(min(mat)));
% muD1_new = normalize(muD1);
muD1_new = muD;
figure(2)
clf
gridsize = 5;
l = (reduceTo - patchsize + 1)^2;
subplot(gridsize, gridsize, 1)
Y_approx = D*(S.*B) + repmat(bias, 1, c.N);
imshow(patch2im(Y_approx(:,1:l), patchsize))
sb = 2;
[~, list_of_f] = sort(sum(B,2));
list_of_f = list_of_f(end:-1:(end - 11));
for i = 1:length(list_of_f)
    sb = i*2;
%     active_f = fs(sb - 1);
    active_f = list_of_f(i);
    tempB = B;
    tempB([1:(active_f - 1), (active_f+1):K1], :) = 0;
    Y_approx = D*(S.*tempB) + repmat(bias, 1, c.N);
    
    subplot(gridsize, gridsize, sb)
    recon = normalize(patch2im(Y_approx(:,1:l), patchsize));
    imshow(recon)
    title(sprintf('Feature %d', active_f))
    
    if totalImages > 1
    
        subplot(gridsize, gridsize, sb + 1)
        recon = normalize(patch2im(Y_approx(:,(l + 1):2*l), patchsize));
        imshow(recon)
        title(sprintf('Feature %d', active_f))
    end
    
end
%% Plot higher level features
normalize = @(mat) (mat - min(min(mat)))/(max(max(mat)) - min(min(mat)));
% muD1_new = normalize(muD1);
figure(2)
clf
gridsize = 5;
l = (reduceTo - patchsize + 1)^2;
subplot(gridsize, gridsize, 1)
Y_approx = D*(S.*B) + repmat(bias, 1, c.N);
imshow(patch2im(Y_approx(:,1:l), patchsize))
sb = 2;
[~, list_of_f] = sort(sum(B2,2));
list_of_f = list_of_f(end:-1:(end - 24));
for i = 1:length(list_of_f)
    sb = i + 1;
%     active_f = fs(sb - 1);
    active_f = list_of_f(i);
    Y2_approx = D2(:, active_f)*(S2(active_f, :).*B2(active_f,:)) + repmat(bias2, 1, c2.N);
    Y2_approx = un_repatch(1./(1+exp(-Y2_approx)), reduceTo, patchsize, K1, totalImages);
%     tempB = Y2_approx > 0.6;
    tempB = binornd(ones(K1, c.N), Y2_approx);
    Y_approx = D*(S.*(tempB.*B)) + repmat(bias, 1, c.N);
    subplot(gridsize, gridsize, sb)
    recon = normalize(patch2im(Y_approx(:,1:l), patchsize));
    imshow(recon)
    title(sprintf('Feature %d', active_f))
%     if totalImages >1
%         subplot(gridsize, gridsize, sb + 1)
%         recon = normalize(patch2im(Y_approx(:,(l + 1):2*l), patchsize));
%         imshow(recon)
%         title(sprintf('Feature %d', active_f))
%     end
    
end
%% Draw Lower Level Patches
figure(5)
clf
% reshape(D2(:, 1), patchsize, patchsize)
[~, list_of_f] = sort(sum(B,2));
list_of_f = list_of_f(end:-1:(end - 24));
for i = 1:25
    subplot(5,5,i)
    j = list_of_f(i);
    imshow(reshape(D(:, j), patchsize, patchsize))
    imagesc(reshape(D(:, j), patchsize, patchsize)); colormap;
    title(sprintf('Feature %d', j))
end
%%
figure(5)
clf
% reshape(D2(:, 1), patchsize, patchsize)
[~, list_of_f] = sort(sum(B2,2));
list_of_f = list_of_f(end:-1:(end - 19));
for i = 1:25
    subplot(5,5,i)
    j = list_of_f(i);
    d2j = D2(:, j);
    d2j = 1./(1 + exp(-d2j));
    d2j = d2j > 0.5;
%     d2j'
    da = reshape(D*d2j(1:K1), patchsize, patchsize);
    db = reshape(D*d2j((K1 + 1):2*K1), patchsize, patchsize);
    dc = reshape(D*d2j((2*K1 + 1):3*K1), patchsize, patchsize);
    dd = reshape(D*d2j((3*K1 + 1):4*K1), patchsize, patchsize);
    I = [da, dc;...
        db, dd];
    imshow(I);
    imagesc(I); colormap;
    title(sprintf('Feature %d', j))
end
%%
Y2_approx = D2(:, 17)*(S2(17, :).*B2(17,:)) + repmat(bias2, 1, c2.N);
% Y2_approx = D2*(S2.*B2) + repmat(bias2, 1, c2.N);
figure(7), clf, imagesc(Y2_approx); colorbar
B_approx = un_repatch(1./(1+exp(-Y2_approx)), reduceTo, patchsize, K1, totalImages);
figure(5), clf
subplot(1,2,1)
imagesc(B_approx); colorbar; title('Recon')
subplot(1,2,2)
imagesc(1./(1+exp(-muB_norm))); colorbar; title('Actual')