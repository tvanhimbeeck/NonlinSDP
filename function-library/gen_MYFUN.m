%% f = gen_MYFUN(PARAM1)
% Template for defingin new concave or convex functions objects compatible with the solver
%   - PARAM1 are optional parameters defining the function

function f = gen_MYFUN(PARAM1)
    
    % optional: import/precompute parameters PARAM2 based on PARAM1
    
    f.fun = @(X)   ( MYFUN_fun ( X,PARAM2 ) );
    f.diff = @(X)  ( MYFUN_diff( X,PARAM2 ) );
    f.hess = @(X,V)( MYFUN_hess( X,V,PARAM2 ) ); %(not used in frank-wolfe)
    f.conv = %'convex' or 'concave'
    f.beta = BETA; % beta value (not used in frank-wolfe)
end

function fval = MYFUN_fun(X,PARAM2)
    % Script evaluating function at rho
end

function grad = MYFUN_diff(X,PARAM2)
    % Script evaluating gradient at rho
end

function HX = MYFUN_hess(X,V,PARAM2)
    % Scipt evaluating derivative of the gradient in the direction V
end