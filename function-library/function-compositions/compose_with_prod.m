%% f3 = compose_with_prod( f1,a )
% constructs the matrix function 
%   	f3(X) := a*f1(X)

function f3 = compose_with_prod( f1,a )
    
    f3.fun  = @(X)    ( a*f1.fun( X )   );
    f3.diff = @(X)    ( a*f1.diff( X )  );
    
    if strcmp(f1.input,'matrix')
        f3.hess = @(X,V)  ( a*f1.hess( X,V ));
    elseif strcmp(f1.input,'vector')||strcmp(f1.input,'scalar')
        f3.hess = @(X)  ( a*f1.hess( X ));
    end
    
    f3.beta = f1.beta;
    
    if strcmp(f1.conv,'convex')&& a>=0
        f3.conv = 'convex';
    elseif strcmp(f1.conv,'convex')&& a<0
        f3.conv = 'concave';
    elseif strcmp(f1.conv,'concave')&& a>=0
        f3.conv = 'concave';
    elseif strcmp(f1.conv,'concave')&& a<0
        f3.conv = 'convex';
    end
    
    f3.input = f1.input;

end