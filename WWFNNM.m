% function X = WWFNNM( Y, NSig, m, Par )
function X=WWFNNM( Y,M,sizeX, m, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   The algorithm solves the following optimization problem
%      1/2||Y-UV'||_F_2 + lambda/2(||U||_F_2 + ||V'||_F_2)
warning('off');
Par.rank=4;
%% initialization
[row, col] = size(Y);
U = randn(row, Par.rank);
tempU = zeros(size(U));
VT = randn(Par.rank, col);
tempVT = zeros(size(VT));
maxiter = 6;
OMG = 1.75;
epsl = 1e-5;

seta = mean(mean(Y.^2, 1)); % may be useless
% lambda = Par.lambdac * Par.c * NSig^2 / sqrt(seta) + epsl; %originally
NSig=10;
W = 1 ./ NSig;
% lambda1 = lambda * col / Par.rank;
% lambda2 = lambda * row / Par.rank;

%% Alternate Direction Optimization
f_curr = 0;
maxiterAqr=3
for i = 1:maxiterAqr
    f_prev = f_curr;
    
    % Fix U and update V
    for j = 1:col
        VT(:, j) = (U' * U + (lambda / W(j)^2) * eye(Par.rank)) \ (U' * Y(:, j));
    end
    tempVT = OMG * VT + (1 - OMG) * tempVT;
    VT = tempVT;
    % Fix V and update U
    U = (Y * diag(W.^2) * VT') / (VT * diag(W.^2) * VT' + lambda * eye(Par.rank));
    tempU = OMG * U + (1 - OMG) * tempU;
    U = tempU;
    % energy function
    DT = norm((Y - U * VT) * diag(W), 'fro') ^ 2;
    %     DT = DT(:)'*DT(:);
    RT = lambda * norm(U, 'fro') ^ 2 + lambda * norm(VT, 'fro') ^ 2;
    f_curr = 0.5 * (DT + RT);
    %     fprintf('FNNM Energy, %d th: %2.8f\n', i, f_curr);
    if (abs(f_prev - f_curr) / f_curr < 0.001)
        break;
    end
    
end

%% result
X = U * VT;
X=M(X,2);
X=X + m;
X=reshape(X,sizeX);
end



