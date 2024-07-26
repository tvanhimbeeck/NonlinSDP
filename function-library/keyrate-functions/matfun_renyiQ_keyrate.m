%% f = matfun_renyiQ_keyrate( beta,input,type )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the function
%       h(X) = - sum_{pk} Tr[ (K_kp X^(1/alpha) K_kp')^alpha]
% with rho_p = K_p rho K_p'

function f = matfun_renyiQ_keyrate( beta,krausP,krausP_K )

    P = length(krausP);
    for p = 1:P
        funp{p} = matfun_K( beta,krausP_K{p} );
    end
    f = matfun_lincomb( funp,ones(1,P),krausP );
end
