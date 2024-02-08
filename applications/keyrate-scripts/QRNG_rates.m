%% generate a random D dimensional QRNG 
d = 2;
k = 2;
p = 2;
z = 10;
gamma = 0.1;
beta = 0.0001;

povmtest = rand_povm( 1,z,d );
povmkey = rand_povm( k,p,d );
rho0 = rand_state( d );
q0 = eval_povm( rho0,povmtest );

%% generate BB84 measurements
d = 2;
z = 10;

KrausP = { eye(d)/sqrt(2),eye(d)/sqrt(2) };
KrausKP = { [1 0],[1 1]/sqrt(2);[0 1], [1 -1]/sqrt(2)};

rho0 = rand_state(2);
povmtest = rand_povm( 1,z,d );
q0 = eval_povm( rho0,povmtest );

%%
options.solver = 'interior-point';
options.epsilon = 1e-4;

%% asymptotic key rates
f = matfun_cond_entr_keyrate( {KrausKP,KrausP},'kraus' );
%f = matfun_cond_entr_keyrate( povmkey,'povm' );
%[rho1,fval] = nlsdp_solve( rho0,f,{eye(d)},[1],[],[],options )


%% finite-size key rates

f1 = matfun_rel_entr_keyrate( q0,povmtest );
fQ = compose_with_prod( fQ,-1 );
fH = matfun_renyiH_keyrate( beta,{KrausKP,KrausP},'kraus' );
fQ = matfun_renyiQ_keyrate( beta,{KrausKP,KrausP},'kraus' );
fQ = compose_with_prod( fQ,-1 );
fsum = compose_with_sum( f1,fH,gamma/beta,(1-gamma) );

%% solving
[rho1,fval1] = nlsdp_solve( rho0,fH,povmtest,q0,[],[],options )
%[rho2,fval2] = nlsdp_solve( rho0,fsum,{eye(d)},[1],{},[],options )


