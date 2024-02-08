%% f = matfun_cond_entr_keyrate( input, type )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)

function f = matfun_renyiH_keyrate( beta,input,type )
    g = matfun_renyiQ_keyrate( beta,input,type );
    f = compose_with_log(g,-1/beta);
end