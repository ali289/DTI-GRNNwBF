% %
% function  [X] =  WWNNM( Y,M, m,sizeX )
% NSig = 40;
% Par.mu=1e-6;
% Par.ReWeiIter=3;
% C=0.1;
% % initialize X as Y or ?
% X = Y;
% for i = 1:Par.ReWeiIter
%     Ystar = X + 1/Par.mu * (Y - X) * diag(NSig).^2;
%     [U, SigmaY, V] =   svd(full(Ystar), 'econ'); % Ystar
%     PatNum        = size(Y,2);
%     TempC  = C * sqrt(PatNum) * 2 * NSig(1)^2; % 2/Par.mu * 
%     [SigmaX, svp] = ClosedWNNM(SigmaY, TempC, eps);
%     X =  U(:,1:svp) * diag(SigmaX) * V(:,1:svp)';
% end
%     X=M(X,2);
%     X=X + m;
%     X=reshape(X,sizeX);
% return;

function  [X]  =  WWNNM( Y, M, m, sizeX )
Iter=3;
C=0.1;
NSig = 40;
[U,SigmaY,V] =   svd(full(Y), 'econ');
PatNum        = size(Y, 2);
Temp            =   sqrt(max( diag(SigmaY) .^ 2 - PatNum * NSig(1) ^ 2, 0 ));
W_Sam = C * sqrt(PatNum) * NSig(1) ^ 2;
for i=1:Iter
    %         W_Vec    =   (C * sqrt(PatNum) * NSig .^ 2) ./ ( Temp + eps );               % Weight vector
    W_Vec  = W_Sam ./ ( Temp + eps );
    SigmaX = soft(SigmaY, diag(W_Vec));
    Temp    = diag(SigmaX);
end
X =  U * SigmaX * V' + m;
X=M(X,2);
X=X + m;
X=reshape(X,sizeX);
return;
