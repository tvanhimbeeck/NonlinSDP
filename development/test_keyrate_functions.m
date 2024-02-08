% %% generate a random D dimensional QRNG 
% d = 2;
% k = 2;
% p = 2;
% z = 10;
% gamma = 0.1;
% beta = 0.0001;
% 
% povmtest = rand_povm( 1,z,d );
% povmkey = rand_povm( k,p,d );
% rho0 = rand_state( d );
% q0 = eval_povm( rho0,povmtest );

%% generate BB84 measurements
d = 2;
z = 10;

KrausP = {eye(d)};
KrausKP = { [1 0];[0 1] };
%KrausP = { eye(d)/sqrt(2),eye(d)/sqrt(2) };
%KrausKP = { [1 0],[1 1]/sqrt(2);[0 1], [1 -1]/sqrt(2)};

rho0 = rand_state(2);
%rho0 = eye(2);

beta = 0.01;
f = matfun_cond_entr_keyrate( {KrausKP,KrausP},'kraus' );
%f = matfun_cond_entr_keyrate( povmkey,'povm' );

fH = matfun_renyiH_keyrate( beta,{KrausKP,KrausP},'kraus' );

fQ = matfun_renyiQ_keyrate( beta,{KrausKP,KrausP},'kraus' );

[f.fun(rho0),fH.fun(rho0),fQ.fun(rho0),(1-fQ.fun(rho0))/beta]

betarange = linspace(0.01,4,100);
for i = 1:length(betarange)
    hrange(i) = H(betarange(i),KrausKP,KrausP,rho0)
end
plot(betarange,hrange)
hold on
plot(0,f.fun(rho0),'*')
function h = H(beta,KrausKP,KrausP,rho)
    f = matfun_renyiQ_keyrate( beta,{KrausKP,KrausP},'kraus' );
    h = (1-f.fun(rho))/beta;
end