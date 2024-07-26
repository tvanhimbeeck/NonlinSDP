%% f = matfun_linear(A)
% Linear function f(X) = Tr[X*A]

function f = matfun_linear(A)
        
    f.fun = @(X)   ( linear_fun ( X,A ) );
    f.diff = @(X)  ( linear_diff( X,A ) );
    f.hess = @(X,V)( linear_hess( X,V,A ) );
    f.conv = 'convex';
    f.beta = 0;
end

function fval = linear_fun(X,A)
    fval = real(trace(X*A));
end

function grad = linear_diff(X,A)
    grad = (A+A')/2;
end

function HX = linear_hess(X,V,A)
    HX = 0;
end