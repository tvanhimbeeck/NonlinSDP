function f = gen_cond_entropy_povm( KrausKP,KrausP )
    

    %%%%not finished !!! use kraus version instead
    
    
    %[KrausKP,KrausP,povmP] = gen_kraus( POVM );
    D = length( KrausKP{1,1} );
    G = gen_GellMann( D );
    
    f{1} = @(x)( cond_entropy(      x,KrausKP,KrausP,G ) );
    f{2} = @(x)( cond_entropy_diff( x,KrausKP,KrausP,G ) );
    f{3} = @(x)( cond_entropy_hess( x,KrausKP,KrausP,G ) );
end

function fval = cond_entropy( x,povmKP,povmP,G )
    
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

function g = cond_entropy_diff( x,povmKP,povmP,G )
    
    [K,P] = size(KrausKP);
    delta = 0;
    rho = vec2mat( x,G );
    S = chol(rho);
    for p = 1:P
        delta = delta - (S\tlogt_m( S*povmP{p}*S' ))/S' ;
        for k = 1:K
            delta = delta + (S\tlogt_m( S*povmKP{k,p}*S' ))/S';
        end
    end
    g = mat2vec( delta,G );
end

function Y = tlogt_m(X)
    % t -> t*log(t) for SDP matrices
    [U,D] = eig(X);
    d = max(diag(D),0);
    Y = U* diag( d.*log(d+(0==d)) )*U';
end

function H = cond_entropy_hess( x,KrausKP,KrausP,G )
    [K,P] = size(KrausKP);
    d = length(x);
    H = zeros( d );
    for i = 1:d, for j = 1:i
        for p = 1:P
            H(i,j) = H(i,j)...
                     - trace( KrausP{p}*G{i}*KrausP{p}'*...
                                logm_frechet_pade(...
                                    KrausP{p}*mat2vec(rho)*KrausP{p}',...
                                    KrausP{p}*G{j}*KrausP{p}')...
                             );
            for k = 1:K
                    H(i,j) = H(i,j)...
                        + trace( KrausKP{k,p}*G{i}*KrausKP{k,p}'*...
                                    logm_frechet_pade(...
                                        KrausKP{k,p}*mat2vec(rho)*KrausKP{k,p}',...
                                        KrausKP{k,p}*G{j}*KrausKP{k,p}' )...
                             );
            end
        end
        H(j,i) = H(i,j);
        end
    end
        
        
end
