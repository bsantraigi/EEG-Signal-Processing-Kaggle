%% Create multivariate basis signals
clear
close all
t = 1:1000;
x = 0:1/80:1;
x = x(2:81);
gam = 1./((x).*(1-x));
b1 = zeros(1000, 80);
for i =1:1000
    b1(i, :) = gam*sin(i/100);
end
% mesh(b1)


gam = sin(x*10);
b2 = zeros(1000, 80);
for i =1:1000
    b2(i, :) = gam*square(i/100);
end
% mesh(b2)


x = 0:1/80:1;
x = x(1:80);
gam = square(x*10);
b3 = zeros(1000, 80);
for i =1:1000
    b3(i, :) = gam*square(i/100);
end
% mesh(b3)

%%
X = zeros(1000, 80, 10);
X(:,:, 1) = 1.0*b1 + 0.0*b2 + 0.0*b3;
X(:,:, 1) = 0.3*b1 + 0.0*b2 + 0.4*b3;
X(:,:, 1) = 0.0*b1 + 0.5*b2 + 0.5*b3;
X(:,:, 1) = 0.7*b1 + 0.0*b2 + 0.9*b3;
X(:,:, 1) = 0.0*b1 + 1.0*b2 + 0.0*b3;
X(:,:, 1) = 0.0*b1 + 0.0*b2 + 1.0*b3;
X(:,:, 1) = 0.0*b1 + 0.0*b2 + 0.0*b3;
X(:,:, 1) = 0.5*b1 + 0.5*b2 + 0.5*b3;
X(:,:, 1) = 0.0*b1 + 0.2*b2 + 0.3*b3;
X(:,:, 1) = 0.7*b1 + 0.1*b2 + 0.7*b3;
%%
dictSize = 5;

D = rand(size(X, 1), size(X,2), dictSize);
norm2 = sqrt(sum(sum(D.*D)));
for n = 1:dictSize
    D(:,:,n) = norm2(:,:,n);
end

mag = @(v) sqrt(sum(v.^2, 1));

tnz = 5;
reduceSize = false;

%% K - SVD original
for it = 1:iter
    W = sparseapprox(X, D, 'OMP', 'tnz',tnz);
    X_a = D*W;
    R = X - X_a;
    if(reduceSize)
        useless = find(sum(W~=0, 2) == 0);
        D(:, useless) = [];
        W(useless, :) = [];
        dictSize = size(D, 2);
        disp(['Reduced dictionary size ' , num2str(dictSize)])
    end
%     clf;
%     plot(X(3000:4000,181),'g')
%     hold on
%     plot(X_a(3000:4000, 181),'r')
%     drawnow;
    
    for k=1:dictSize
        I = find(W(k,:));
        if size(I, 2) == 0
            continue
        end
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        % U is normalized
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end
%     er = mean(mag(R));
%     disp(['Error ', num2str(it), ':', num2str(er)])
%     errors(it) = mean(mag(R));
    snrs = snr_2(X, D*W);
    fprintf('Range of SNR %d: (%f, %f)\n', it, min(snrs), max(snrs));
    
%     plot(t, errors);
%     drawnow;
    
end
W = sparseapprox(X, D, 'OMP', 'tnz',tnz);


