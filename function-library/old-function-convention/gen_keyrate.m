%% f = gen_keyrate( input, type )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Builds the keyrate function object, which can be specified in terms of POVM elements 
%   - type = 'povm'
%   - input = two index cell containing POVM elements
% or in terms of Kraus operators
%   - type = 'kraus'
%   - input = two index cell containing Kraus operators
% For non-full rank POVM's Kraus operators
% can be computed but for non full rank POVM
%
function f = gen_keyrate( input, type )
    
    if strcmp(type,'povm') %POVM type input
        [KrausKP,KrausP] = gen_kraus( input );
    elseif strcmp(type,'kraus') % KRAUS type input
        KrausKP = input{1};
        KrausP = input{2};
    end
    
    load tau_r8_zeros
    
    f.fun = @(rho)   ( cond_entropy_fun(      rho,  KrausKP,KrausP ) );
    f.diff = @(rho)  ( cond_entropy_diff( rho,  KrausKP,KrausP ) );
    f.hess = @(rho,X)( cond_entropy_hess( rho,X,KrausKP,KrausP,pzero,qzero ) );
    f.conv = 'convex';
    f.beta = 1;
end

function fval = cond_entropy_fun( rho,KrausKP,KrausP )
    
    [K,P] = size(KrausKP);
    fval = 0;
    for p = 1:P
        fval = fval - VNent( KrausP{p}*rho*KrausP{p}' );
        for k = 1:K
            fval = fval + VNent( KrausKP{k,p}*rho*KrausKP{k,p}' );
        end
    end
    fval = real(fval);
end

function grad = cond_entropy_diff( rho,KrausKP,KrausP )
    
    [K,P] = size(KrausKP);
    grad = 0;
    for p = 1:P
        grad = grad + KrausP{p}'*logm( KrausP{p}*rho*KrausP{p}' )*KrausP{p};
        for k = 1:K
            grad = grad - KrausKP{k,p}'*logm( KrausKP{k,p}*rho*KrausKP{k,p}' )*KrausKP{k,p};
        end
    end
    grad = (grad+grad')/2;
    
end

function HX = cond_entropy_hess( rho,X,KrausKP,KrausP,pzero,qzero )
    [K,P] = size(KrausKP);
    D = length(rho);
    HX = zeros( D );
    for p = 1:P
        HX = HX + KrausP{p}'*...
                        logm_frechet_pade_herm(...
                            KrausP{p}*rho*KrausP{p}',...
                            KrausP{p}*X*KrausP{p}',...
                            pzero,qzero)*...
                        KrausP{p};
        for k = 1:K
            HX = HX - KrausKP{k,p}'*...
                            logm_frechet_pade_herm(...
                                KrausKP{k,p}*rho*KrausKP{k,p}',...
                                KrausKP{k,p}*X*KrausKP{k,p}',...
                                pzero,qzero)*...
                            KrausKP{k,p};
        end
    end
    HX = (HX + HX')/2;
end
