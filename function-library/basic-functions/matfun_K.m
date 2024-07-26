%% f = matfun_K( beta,K )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the concave function
%       h(X) = sum_{x} Tr[ (K_x rho^(1/alpha) K_x')^alpha]
% with K a family of kraus operators

function f = matfun_K( beta,K )

    f.fun = @(rho)   ( fun ( beta,K,rho ) );
    f.diff = @(rho)  ( diff( beta,K,rho ) );
    f.hess = @(rho,X)( hess( beta,K,rho,X ) );
    f.conv = 'concave';
    f.beta = 1;
    f.input = 'matrix';
end

function fval = fun( beta,K,rho )
    
    % sum_x Tr [(K rho^(1/alpha) K' )^alpha ]
    fval = 0;
    alpha = 1+beta;
    rho_beta = mpower(rho,1/alpha);
    for k = 1:length(K)
        fval = fval + trace( mpower( K{k}*rho_beta*K{k}',alpha ) );
    end
    fval = real(fval);
end

function grad = diff( beta,K,rho )
    % sum d()^1/alpha|_rho[ sum K(K rho^1/alpha K)^beta K ]
    alpha = 1+beta;
    rho_beta = mpower(rho,1/alpha);
    Delta = 0;
    for k = 1:length(K)
        Delta = Delta + K{k}'*mpower(K{k}*rho_beta*K{k}',beta)*K{k};
    end
    grad = alpha*power_grad(1/alpha,rho,Delta);
    grad = (grad+grad')/2;
end

function HX = hess( beta,K,rho,X )
    % Delta = sum K(K rho^1/alpha K)^beta K
    % d2()^1/alpha|_rho[X,Delta]
    % + sum d()^1/alpha|_rho[ K d()^beta|_(K rho^1/alpha K) [ K d()^1/alpha|_rho[X] K] K] 
    alpha = 1+beta;
    X = (X+X')/2;
    rho_beta = mpower( rho,1/(1+beta) );
    Delta = 0;
    for k = 1:length(K)
        Delta = Delta + K{k}'*mpower( K{k}*rho_beta*K{k}',beta )*K{k};
    end
    HX = power_hess( 1/alpha,rho,X,Delta );
    
    V = power_grad( 1/alpha,rho,X );
    for k = 1:length(K)
        temp = K{k}*V*K{k}';
        temp2 = K{k}'*power_grad( beta,K{k}*rho_beta*K{k}',temp )*K{k};
        HX = HX + power_grad( 1/alpha,rho,temp2 );
    end
    HX = alpha*(HX + HX')/2;
end