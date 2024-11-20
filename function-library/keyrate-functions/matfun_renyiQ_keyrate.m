%% f = matfun_renyiQ_keyrate( beta,krausP,krausP_K )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the function
%       h(rho) = sum_{pk} Tr[ (krausP_K{k}{p} rho_p^(1/alpha) krausP_K{k}{p}')^alpha]
% with rho_p = krausP{p} rho krausP{p}' and alpha = 1+beta

function f = matfun_renyiQ_keyrate( beta,krausP,krausP_K )

    P = length(krausP);
    for p = 1:P
        funp{p} = matfun_K( beta,krausP_K{p} );
    end
    f = compose_lincomb( funp,ones(1,P),krausP );
end