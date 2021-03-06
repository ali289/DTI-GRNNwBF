function [X]  = SVS(y,M,sizeX,rankX, lambda)


% Copyright (c) Angshul Majumdar 2010

    err = 1e-12;
    x_initial = zeros(prod(sizeX),1); % leave it as it is
    normfac = 1; % leave it as it is
    insweep = 50;%20;
    tol = 1e-4;%1e-4;    
    decfac = 0.7;%earlier it was 0.9

alpha =1.1*normfac;
x = x_initial;
% lambdaInit = decfac*max(abs(M(y,2))); lambda = lambdaInit;
f_current = norm(y-M(x,1)) + lambda*norm(x,1);

xs=[]; xs2=[];
% while lambda >lambdaInit*tol
    %lambda
    for ins = 1:insweep
        
        f_previous = f_current;
        b = x + (1/alpha)*M(y - M(x,1),2); % Landweber Iteration
        B = reshape(b,sizeX);
        try
        [U,S,V] = svd(B,'econ'); %lansvd(B,rankX,'L'); %  uses more efficient lansvd from PROPACK
        
        s = SoftTh(diag(S),lambda/(2*alpha));
        S = diag(s);
        X = U*S*V';
        X(X<0)=0;
        x = X(:);
        f_current = norm(y-M(x,1)) + lambda*norm(x,1);
        
         xs=[xs; norm(y-M(x,1),'fro')];
        xs2=[xs2; norm(x,1)];
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
        
        catch
            fprintf('SVD didnt converge!')
        end
        
    end
    
    if norm(y-M(x,1))<err
      %  break;
    end
   
    lambda = decfac*lambda;
    
    
         
% end

    function  z = SoftTh(s,thld)
        z = sign(s).*max(0,abs(s)-thld); 
    end
end
