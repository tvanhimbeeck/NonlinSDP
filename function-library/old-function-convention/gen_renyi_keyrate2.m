%% f = gen_renyi_keyrate2( beta,POVM,Op,gamma )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% f(rho) = (1-gamma) sum_{ax} Tr[ (K_ax rho_x^(1/alpha) K_kp')^alpha]
%          + gamma   sum_z p(z) Tr[Op rho]
% with rho_x = K_x rho K_x'

function f = gen_renyi_keyrate2( beta,POVM,Op,gamma )

    [KrausKP,KrausP] = gen_kraus( POVM );
        
    f.fun  = @(rho)  ( renyi_keyrate2_fun ( beta, rho,  KrausKP,KrausP,Op,gamma ) );
    f.diff = @(rho)  ( renyi_keyrate2_diff( beta, rho,  KrausKP,KrausP,Op,gamma ) );
    f.hess = @(rho,X)( renyi_keyrate2_hess( beta, rho,X,KrausKP,KrausP,Op,gamma ) );
    f.conv = 'concave';
    f.beta = 1;
end

function fval = renyi_keyrate2_fun( beta,rho,KrausKP,KrausP,Op,gamma )
    
    [K,P] = size(KrausKP);
    fval = 0;
    for p = 1:P
        rho_p = KrausP{p}*rho*KrausP{p}';
        rho_p_beta = mpower(rho_p,1/(1+beta));
        for k = 1:K
            fval = fval + trace( mpower(KrausKP{k,p}*rho_p_beta*KrausKP{k,p}',1+beta) );
        end
    end
    fval = (1-gamma)*fval;
    fval = fval + gamma*trace(rho*Op);
    fval = real(fval);
end

function grad = renyi_keyrate2_diff( beta,rho,KrausKP,KrausP,Op,gamma )
    
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
    grad = (1-gamma)*grad + gamma*Op;
    grad = (grad+grad')/2;
end

function HX = renyi_keyrate2_hess( beta,rho,X,KrausKP,KrausP,Op,gamma )
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
    HX = (1-gamma)*alpha*(HX + HX')/2;
end