function dfval = cond_entropy_diff( rho,keyMap )
    % note : base = 2
    
    D = length(rho);
    [K,P] = size( keyMap );
    [KrausKP,KrausP] = gen_kraus( keyMap );
    
    dfval = 0;
    for p = 1:P
        dfval = dfval + KrausP{p}*my_logm( KrausP{p}*rho*KrausP{p} )*KrausP{p};
        for k = 1:K
            dfval = dfval - KrausKP{k,p}*my_logm( KrausKP{k,p}*rho*KrausKP{k,p} )*KrausKP{k,p};
        end
    end
    dfval = 1/2*(dfval + dfval')/log(2);
end

function logM = my_logm( M )
    M = (M + M')/2;
    logM = logm( (1 - 1e-10)*M + 1e-10*eye( length( M ) ));
end