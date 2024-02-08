function H = cond_entropy_hess( rho,keyMap,G )

    % returns the hessian of cond_entropy as a matrix 
    % in the orthonormal basis G;
    
    [K,P] = size( keyMap );
    [KrausP,KrausKP] = gen_kraus( keyMap );
    
    H = 0;
    for p = 1:P
        K = KrausP{p};
        H = H - VNent_hess( K*rho*K,G,K );
        for k = 1:K
            K = KrausKP{k,p};
            H = H + VNent_hess( K*rho*K,G,K );
        end
    end
end