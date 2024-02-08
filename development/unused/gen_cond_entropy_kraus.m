function f = gen_cond_entropy_kraus( KrausKP,KrausP )

    D = length( KrausKP{1,1} );
    G = gen_GellMann( D );
    load tau_r8_zeros
    f.fun = @(x)( cond_entropy(      x,KrausKP,KrausP,G ) );
    f.diff = @(x)( cond_entropy_diff( x,KrausKP,KrausP,G ) );
    f.hess = @(x)( cond_entropy_hess( x,KrausKP,KrausP,G,pzero,qzero ) );
    f.beta = 1;
end

function fval = cond_entropy( x,KrausKP,KrausP,G )
    
    [K,P] = size(KrausKP);
    fval = 0;
    rho = vec2mat( x,G );
    for p = 1:P
        fval = fval - VNent( KrausP{p}*rho*KrausP{p}' );
        for k = 1:K
            fval = fval + VNent( KrausKP{k,p}*rho*KrausKP{k,p}' );
        end
    end
    fval = real(fval);
end

function g = cond_entropy_diff( x,KrausKP,KrausP,G )
    
    [K,P] = size(KrausKP);
    delta = 0;
    rho = vec2mat( x,G );
    for p = 1:P
        delta = delta + KrausP{p}'*logm( KrausP{p}*rho*KrausP{p}' )*KrausP{p};
        for k = 1:K
            delta = delta - KrausKP{k,p}'*logm( KrausKP{k,p}*rho*KrausKP{k,p}' )*KrausKP{k,p};
        end
    end
    g = mat2vec( delta,G );
end

function H = cond_entropy_hess( x,KrausKP,KrausP,G,pzero,qzero )
    [K,P] = size(KrausKP);
    d = length(x);
    H = zeros( d );
    rho = vec2mat(x,G);
    for i = 1:d
        HGi = 0;
        for p = 1:P
            HGi = HGi + KrausP{p}'*...
                            logm_frechet_pade_herm(...
                                KrausP{p}*rho*KrausP{p}',...
                                KrausP{p}*G{i}*KrausP{p}',...
                                pzero,qzero)*...
                            KrausP{p};
            for k = 1:K
                HGi = HGi - KrausKP{k,p}'*...
                                logm_frechet_pade_herm(...
                                    KrausKP{k,p}*rho*KrausKP{k,p}',...
                                    KrausKP{k,p}*G{i}*KrausKP{k,p}',...
                                    pzero,qzero)*...
                                KrausKP{k,p};
                for j = 1:d
                    H(j,i) = real(inner_prod( G{j},HGi ));
                end
            end
        end
    end
    H = (H + H')/2;
end
