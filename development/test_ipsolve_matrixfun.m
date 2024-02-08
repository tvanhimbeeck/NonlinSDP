d = 10;
z = 10;

%KrausP = {eye(d)};
%KrausKP = { [1 0];[0 1] };
K = 10; P = 3;
povmkey = rand_povm(K,P,d);

X0 = rand_state(d);
povmtest = rand_povm( 1,z,d );
q0 = eval_povm( X0,povmtest );
X0 = eye(d)/d;
%X0 = rand_state(d);
%%
beta_range = 0.1;%0.01:0.05:0.99;
n = length(beta_range);
key_range = [];
for i = 1:length(beta_range)
    beta = beta_range(i);
    f = matfun_rel_entr_keyrate( q0,povmtest ); % (d,d) -> 1 convex function

    %F = matfun_renyiQ_keyrate( beta,{KrausKP,KrausP},'kraus' ); 
    F = matfun_renyiQ_keyrate( beta,povmkey,'povm' ); 
    F = compose_with_prod( F,-1 ); % (d,d) -> 1 convex function

    g = fun_log();
    g = compose_with_prod( g,-1/beta );
    g = precompose_with_prod( g,-1 ); % 1 -> 1 convex function

    options.solver = 'interior-point';
    options.epsilon = 1e-5;
    options.verbose = 'real-time';
    [X1,fval,output] = ipsolve_matrixset( X0,f,g,F,{eye(d)},[1],options );
    [f.fun(X0),g.fun(F.fun(X0)),f.fun(X1),g.fun(F.fun(X1))]
    key_range = [key_range, fval];
end
plot(beta_range,key_range)

