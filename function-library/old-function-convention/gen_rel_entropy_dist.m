%% f = gen_rel_entropy_dist(q0,POVM)
% Generates the matrix function
%       f(X) = D( q[X]||q0 )
% where 
%       - D(q||p) is the relative entropy for probability distributions 
%       - q_i[X] = Tr[ POVM{i}*X ]

function f = gen_rel_entropy_dist(q0,POVM)
  
    f.fun = @(X)   ( fun ( X,q0,POVM ) );
    f.diff = @(X)  ( diff( X,q0,POVM ) );
    f.hess = @(X,V)( hess( X,V,POVM ) );
    f.conv = 'convex';
    f.beta = 1;
end

function fval = fun(X,q0,POVM)
    for i = 1:length(POVM)
        q(i) = trace(POVM{i}*X);
    end
    fval = sum(q.*(log(q) - log(q0)));
end

function grad = diff(X,q0,POVM)
    grad = 0;
    for i = 1:length(POVM)
        q(i) = trace(POVM{i}*X);
        grad = grad + (log(q(i))-log(q0(i))+1)*POVM{i};
    end
end

function HX = hess(X,V,POVM)
    HX = zeros(length(X));
    for i = 1:length(POVM)
        q(i) = trace(POVM{i}*X);
        v(i) = trace(POVM{i}*V);
        HX = HX + POVM{i}*v(i)/q(i);
    end
end