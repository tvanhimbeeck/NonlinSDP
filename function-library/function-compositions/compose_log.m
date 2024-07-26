%% f2 = compose_log( f1,a )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)
%
% Constructs the convex/concave matrix function 
%   	f2(X) := a*log(f1(X)) 
% where 
%       - f1 is a concave matrix function
%       - a is a real coefficient

function f2 = compose_log( f1,a )
    
    if ~strcmp( f1.conv,'concave' )
        f2 = [];
        fprintf( 'forbidden composition\n' )
    elseif strcmp( f1.conv,'concave' )
        f2.fun = @(X)( a*log( f1.fun( X ) ) );
        f2.diff = @(X)( a*f1.diff(X)/f1.fun(X) );
        f2.hess = @(X,V) ( a*( f1.hess(X,V)/f1.fun(X) - f1.diff(X)*inner_prod( V,f1.diff(X) )/(f1.fun(X)^2 ) ) );
        f2.beta = 0; %% not concordant
        if a >= 0
            f2.conv = 'concave';
        elseif a < 0
            f2.conv = 'convex';
        end
    end
end