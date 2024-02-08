%% f = gen_trlog()

function f = gen_trlog()
        
    f.fun = @(X)   ( trlog_fun ( X) );
    f.diff = @(X)  ( trlog_diff( X ) );
    f.hess = @(X,V)( trlog_hess( X,V ) );
    f.conv = 'convex';
    f.beta = 2;
end

function fval = trlog_fun( X )
    fval = -trace(log(X));
end

function grad = trlog_diff( X )
    grad = -inv(X);
end

function HX = trlog_hess( X,V )
    HX = inv(X)*V*inv(X);
end