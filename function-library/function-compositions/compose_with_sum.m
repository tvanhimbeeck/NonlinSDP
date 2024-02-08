%% f3 = compose_with_sum( f1,f2,a,b )
% constructs the matrix function 
%   	f3(X) := a*f1(X) + b*f2(X) 
% if the convexities match

function f3 = compose_with_sum( f1,f2,a,b )
    if strcmp(f1.conv,'concave')
        c1 = -1;
    elseif strcmp(f1.conv,'convex')
        c1 = 1;
    end
    if strcmp(f2.conv,'concave')
        c2 = -1;
    elseif strcmp(f2.conv,'convex')
        c2 = 1;
    end
    
    if a*c1*b*c2 < 0 %ie. convexities don't match
        fprintf( 'forbidden composition\n' )
        
    elseif a*c1*b*c2 >= 0
        f3.fun  = @(X)    ( a*f1.fun( X )   + b*f2.fun( X )   );
        f3.diff = @(X)    ( a*f1.diff( X )  + b*f2.diff( X )  );
        f3.hess = @(X,V)  ( a*f1.hess( X,V )+ b*f2.hess( X,V ));
        f3.beta = max( f1.beta,f2.beta );
        if a*c1 >= 0
            f3.conv = 'convex';
        elseif a*c1 < 0
            f3.conv = 'concave';
        end
    end
end