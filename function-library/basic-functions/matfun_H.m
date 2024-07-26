%% f = matfun_H( K )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the concave function
%       f(X) =  sum_{x} S(K_x X K_x') - S(X)
% with K a family of kraus operators such that sum_k K{k}'*K{k} = id
%
function f = matfun_H( K )
    
    load tau_r8_zeros
    
    f.fun = @(rho)   ( fun(  K,rho ) );
    f.diff = @(rho)  ( diff( K,rho ) );
    f.hess = @(rho,X)( hess( K,rho,X,pzero,qzero ) );
    f.conv = 'concave';
    f.beta = 1;
end

function fval = fun( K,rho )
    
    fval = - VNent( rho_p );
    for k = 1:length(K)
        fval = fval + VNent( K{k}*rho*K{k}' );
    end
    fval = real(fval);
end

function grad = diff( K,rho )
    
    grad = logm( rho );
    for k = 1:length(K)
        grad = grad - K{k}'*logm( K{k}*rho*K{k}' )*K{k};
    end
    grad = (grad+grad')/2;
    
end

function HX = hess( K,rho,X,pzero,qzero )
    
    HX = logm_frechet_pade_herm( rho,X,pzero,qzero );
    for k = 1:length(K)
        HX = HX - K{k}'*logm_frechet_pade_herm(K{k}*rho*K{k}',K{k}*X*K{l}',pzero,qzero)*K{k};
    end
    HX = (HX + HX')/2;
end
