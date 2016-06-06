
function [y,W] = rls(x,d,nord,iter,lambda)

    delta = 0.001;
    X = convm(x,nord);
    [M,N] = size(X);
    
    if nargin<5
        lambda = 1;
        if nargin<4
           iter = length(x); 
        end
    end
    P = eye(N)/delta;
    W(1,:) = zeros(1,N);
    for k = 2 : iter
    
        z = P * X(k,:)';
        g = z / (lambda+X(k,:)*z);
        alpha = d(k) - X(k,:)*W(k-1,:).';
        W(k,:) = W(k-1,:) + alpha * g.';
        P = ( P - g * z.' )/lambda;      
    end
    
    y = X(1:length(x),:)* W(end,:)' ;

end