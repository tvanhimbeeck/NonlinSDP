% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function [KrausKP,KrausP] = gen_kraus( povm )
    
    [K,P] = size( povm );
    KrausKP = cell(K,P);
    KrausP = cell(1,P);
    for p = 1:P
        M{p} = 0;
        for k = 1:K
            M{p} = M{p} + povm{ k,p };
        end
        KrausP{p} = chol( M{p} );
    end
    
    for p = 1:P
        for k = 1:K
            KrausKP{k,p} = chol( povm{ k,p } )*inv(KrausP{p});
        end
    end
end