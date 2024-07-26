%% f = matfun_K( beta,K )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the concave function
%       h(X) = sum_{x} Tr[ (K_x rho^(1/alpha) K_x')^alpha]
% with K a family of kraus operators

function f = matfun_K( beta,kraus )

    f.fun = @(rho)   ( fun ( beta, rho,  kraus ) );
    f.diff = @(rho)  ( diff( beta, rho,  kraus ) );
    f.hess = @(rho,X)( hess( beta, rho,X,kraus ) );
    f.conv = 'concave';
    f.beta = 1;
    f.input = 'matrix';
end

function fval = fun( beta,rho,kraus )
    
    % sum_x Tr [(K rho^(1/alpha) K' )^alpha ]
    fval = 0;
    alpha = 1+beta;
    rho_beta = mpower(rho,1/alpha);
    for k = 1:length(kraus)
        fval = fval + trace( mpower( kraus{k}*rho_beta*kraus{k}',alpha ) );
    end
    fval = real(fval);
end

function grad = diff( beta,rho,kraus )
    % sum d()^1/alpha|_rho[ sum K(K rho^1/alpha K)^beta K ]
    alpha = 1+beta;
    rho_beta = mpower(rho,1/alpha);
    Delta = 0;
    for k = 1:length(kraus)
        Delta = Delta + kraus{k}'*mpower(kraus{k}*rho_beta*kraus{k}',beta)*kraus{k};
    end
    grad = alpha*power_grad(1/alpha,rho,Delta);
    grad = (grad+grad')/2;
end

function HX = hess( beta,rho,X,kraus )
    % Delta = sum K(K rho^1/alpha K)^beta K
    % d2()^1/alpha|_rho[X,Delta]
    % + sum d()^1/alpha|_rho[ K d()^beta|_(K rho^1/alpha K) [ K d()^1/alpha|_rho[X] K] K] 
    alpha = 1+beta;
    X = (X+X')/2;
    rho_beta = mpower( rho,1/(1+beta) );
    Delta = 0;
    for k = 1:length(kraus)
        Delta = Delta + kraus{k}'*mpower( kraus{k}*rho_beta*kraus{k}',beta )*kraus{k};
    end
    HX = power_hess( 1/alpha,rho,X,Delta );
    
    V = power_grad( 1/alpha,rho,X );
    for k = 1:length(kraus)
        temp = kraus{k}*V*kraus{k}';
        temp2 = kraus{k}'*power_grad( beta,kraus{k}*rho_beta*kraus{k}',temp )*kraus{k};
        HX = HX + power_grad( 1/alpha,rho,temp2 );
    end
    HX = alpha*(HX + HX')/2;
end