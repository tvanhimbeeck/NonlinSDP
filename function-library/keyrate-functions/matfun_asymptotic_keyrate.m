%% f = matfun_asymptotic_keyrate( krausP,krausP_K )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Computes the concave asymptotic key rate function
%       sum_{kp} S(krausP_K{p}{k}*krausP{p}*rho*krausP{p}'*krausP_K{p}{k}') )
%           - sum_{p} S(krausP{p}*rho*krausP{p})

function f2 = matfun_asymptotic_keyrate( krausP,krausP_K )
    
    P = length(krausP);
    for p = 1:P
        f{p} = matfun_H( krausP_K );
    end
    f2 = compose_lincomb( f,ones(1,P),krausP );
    
end