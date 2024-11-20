%% f = matfun_weightedrenyiQ_keyrate( beta,protocol,f,f0 )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Implements the function
%       h(rho) = gamma sum_z exp(beta f(z)) trace[M_z rho]
%                + (1-gamma) sum_{pk} Tr[ (krausP_K{k}{p} rho_p^(1/alpha) krausP_K{k}{p}')^alpha]
% with rho_p = krausP{p} rho krausP{p}' and alpha = 1+beta

function fun = matfun_weightedrenyiQ_keyrate( beta,protocol,f,f0 )
    
    % test round
    if strcmp( protocol.test,"povm" )
        povmZ = protocol.povmZ;
        X = cellvecinner(povmZ,exp(beta*f));
    elseif strcmp( protocol.test,"observables" )
        Xfun = protocol.Xfun;
        X = Xfun(beta*(1-gamma)*y/gamma);
    end
    fX = matfun_linear(X);
    
    % key round
    krausP = protocol.krausP;
    krausP_K = protocol.krausP_K;
    fQ = matfun_renyiQ_keyrate( beta,krausP,krausP_K );
    
    % linear combination
    gamma = protocol.gamma;
    fun = compose_lincomb( {fX,fQ},[gamma,(1-gamma)*exp(beta*f0)] );
end