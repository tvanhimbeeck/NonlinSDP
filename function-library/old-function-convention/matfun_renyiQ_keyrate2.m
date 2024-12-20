%% f = matfun_renyiQ_keyrate( beta,input,type )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the function
%       h(X) = - sum_{pk} Tr[ (K_kp rho_p^(1/alpha) K_kp')^alpha]
% with rho_p = K_p rho K_p'

function f = matfun_renyiQ_keyrate( beta,input,type )

    if strcmp(type,'povm') %POVM type input
        [KrausKP,KrausP] = gen_kraus( input );
    elseif strcmp(type,'kraus') % KRAUS type input
        KrausKP = input{1};
        KrausP = input{2};
    end
        
    f.fun = @(rho)   ( fun ( beta, rho,  KrausKP,KrausP ) );
    f.diff = @(rho)  ( diff( beta, rho,  KrausKP,KrausP ) );
    f.hess = @(rho,X)( hess( beta, rho,X,KrausKP,KrausP ) );
    f.conv = 'concave';
    f.beta = 1;
    f.input = 'matrix';
end

function fval = fun( beta,rho,KrausKP,KrausP )
    
    [K,P] = size(KrausKP);
    fval = 0;
    for p = 1:P
        rho_p = KrausP{p}*rho*KrausP{p}';
        rho_p_beta = mpower(rho_p,1/(1+beta));
        for k = 1:K
            fval = fval + trace( mpower( KrausKP{k,p}*rho_p_beta*KrausKP{k,p}',1+beta ) );
        end
    end
    fval = real(fval);
end

function grad = diff( beta,rho,KrausKP,KrausP )
    
    [K,P] = size(KrausKP);
    grad = 0;
    for p = 1:P
        rho_p = KrausP{p}*rho*KrausP{p}';
        rho_p_beta = mpower(rho_p,1/(1+beta));
        for k = 1:K
            X = KrausKP{k,p}'*mpower(KrausKP{k,p}*rho_p_beta*KrausKP{k,p}',beta)*KrausKP{k,p};
            grad = grad + (1+beta)*power_grad(1/(1+beta),rho_p,X);
        end
    end
    grad = (grad+grad')/2;
end

function HX = hess( beta,rho,X,KrausKP,KrausP )
    [K,P] = size(KrausKP);
    D = length(rho);
    HX = zeros( D );
    alpha = 1+beta;
    for p = 1:P
        rho_p = KrausP{p}*rho*KrausP{p}';
        X_p = KrausP{p}*X*KrausP{p}';
        A = 0;
        for k = 1:K % sum_k K(K rho^(1/alpha) K)^beta K
            A = A + KrausKP{k,p}'*mpower( KrausKP{k,p}*mpower( rho_p,1/(1+beta) )*KrausKP{k,p}',beta )*KrausKP{k,p};
        end
        A_p = KrausP{p}'*A*KrausP{p};
        HX = HX + power_hess( 1/alpha,rho_p,X_p,A_p );
        V = power_grad( 1/alpha,rho_p,X_p );
        for k = 1:K
            temp = KrausKP{k,p}*V*KrausKP{k,p}';
            temp2 = KrausKP{k,p}'*power_grad( beta,KrausKP{k,p}*mpower(rho_p,1/alpha)*KrausKP{k,p}',temp )*KrausKP{k,p};
            HX = HX + KrausP{p}'*power_grad( 1/alpha,rho_p,temp2 )*KrausP{p};
        end
    end
    HX = alpha*(HX + HX')/2;
end