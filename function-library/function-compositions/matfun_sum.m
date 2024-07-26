%% f3 = matfun_lincomb( f,a )
% constructs the matrix function 
%   	f3(X) := sum_i a(i) * f{i}(X) 
% if the convexities match

function f3 = matfun_lincomb( f,a )
    
    k = lenght(f);
    % check matching convexities
    for i = 1:k   
        if strcmp(f1.conv,'concave')
            c(i) = -a(i);
        elseif strcmp(f1.conv,'convex')
            c(i) = a(i);
        end
    end
    if 
    
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