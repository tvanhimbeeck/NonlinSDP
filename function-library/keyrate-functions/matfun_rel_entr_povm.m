%% f = matfun_rel_entr_keyrate( q0,POVM )
% Generates the matrix function
%       f(X) = D( q0||q[X] )
% where 
%       - D(q||p) is the relative entropy for probability distributions 
%       - q_i[X] = Tr[ POVM{i}*X ]

function f = matfun_rel_entr_povm( q0,POVM )
  
    f.fun = @(X)   ( fun ( X,q0,POVM ) );
    f.diff = @(X)  ( diff( X,q0,POVM ) );
    f.hess = @(X,V)( hess( X,V,POVM ) );
    f.conv = 'convex';
    f.beta = 1;
    f.input = 'matrix'; % not used
end

function fval = fun( X,q0,POVM )
    for i = 1:length(POVM)
        q(i) = trace(POVM{i}*X);
    end
    fval = sum(q0.*(log(q0) - log(q)));
    fval = real(fval);
end

function grad = diff( X,q0,POVM )
    grad = 0;
    for i = 1:length(POVM)
        q(i) = trace(POVM{i}*X);
        grad = grad - q0(i)/q(i)*POVM{i};
    end
    grad = (grad+grad')/2;
end

function HX = hess( X,V,POVM )
%       To be completed

%     HX = zeros(length(X));
%     for i = 1:length(POVM)
%         q(i) = trace(POVM{i}*X);
%         v(i) = trace(POVM{i}*V);
%         HX = HX + POVM{i}*v(i)/q(i);
%     end
    HX = NaN;
end