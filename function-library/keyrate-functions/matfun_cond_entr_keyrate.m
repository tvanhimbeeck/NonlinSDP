%% f = matfun_cond_entr_keyrate( input, type )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Builds the keyrate function object, which can be specified in terms of POVM elements 
%   - type = 'povm'
%   - input = two index cell containing POVM elements
% or in terms of kraus operators
%   - type = 'kraus'
%   - input = two index cell containing kraus operators
% For non-full rank POVM's kraus operators
% can be computed but for non full rank POVM
%
function f = matfun_cond_entr_keyrate( input, type )
    
    if strcmp(type,'povm') %POVM type input
        [krausP_K,krausP] = gen_kraus( input ); %% needs modification
    elseif strcmp(type,'kraus') % kraus type input
        krausP = input{1};
        krausP_K = input{2};
    end
    d = size(krausP,2);
    
    load tau_r8_zeros
    
    f.fun = @(rho)   ( cond_entropy_fun(  rho,  krausP_K,krausP ) );
    f.diff = @(rho)  ( cond_entropy_diff( rho,  krausP_K,krausP ) );
    f.hess = @(rho,X)( cond_entropy_hess( rho,X,krausP_K,krausP,pzero,qzero ) );
    f.conv = 'convex';
    f.beta = 1;
end

function fval = cond_entropy_fun( rho,krausP_K,krausP )
    
    P = length(krausP);
    fval = 0;
    for p = 1:P
        rho_p = krausP{p}*rho*krausP{p}';
        fval = fval - VNent( rho_p );
        for k = 1:length(krausP_K{p})
            fval = fval + VNent( krausP_K{p}{k}*rho_p*krausP_K{p}{k}' );
        end
    end
    fval = real(fval);
end

function grad = cond_entropy_diff( rho,krausP_K,krausP )
    
    P = length(krausP_K);
    grad = 0;
    for p = 1:P
        rho_p = krausP{p}*rho*krausP{p}';
        grad = grad + krausP{p}'*logm( rho_p )*krausP{p};
        for k = 1:length(krausP_K{p})
            grad = grad - krausP{p}'*krausP_K{p}{k}'*logm( krausP_K{p}{k}*rho_p*krausP_K{p}{k}' )*krausP_K{p}{k}*krausP{p};
        end
    end
    grad = (grad+grad')/2;
    
end

function HX = cond_entropy_hess( rho,X,krausP_K,krausP,pzero,qzero )
    
    P = length(krausP_K);
    D = length(rho);
    HX = zeros( D );
    for p = 1:P
        HX = HX + krausP{p}'*...
                        logm_frechet_pade_herm(...
                            krausP{p}*rho*krausP{p}',...
                            krausP{p}*X*krausP{p}',...
                            pzero,qzero)*...
                        krausP{p};
        for k = 1:length(krausP_K{p})
            HX = HX - krausP{p}'*krausP_K{p}{k}'*...
                            logm_frechet_pade_herm(...
                                krausP_K{p}{k}*krausP{p}*rho*krausP{p}'*krausP_K{p}{k}',...
                                krausP_K{p}{k}*krausP{p}*X*krausP{p}'*krausP_K{p}{k}',...
                                pzero,qzero)*...
                            krausP_K{p}{k}*krausP{p};
        end
    end
    HX = (HX + HX')/2;
end
