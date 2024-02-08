% test of rel entropy function
K = 20;
POVM = rand_keyMap(K,1,4);
S0 = rand_state(4);
S1 = rand_state(4);
for i = 1:K
    q0(i) = trace(S0*POVM{i});
end
f = gen_rel_entropy_dist(q0,POVM);

options.solver = 'frank-wolfe';
options.epsilon = 1e-3;
[S2,obj] = nlsdp_solve( S1,f,{eye(4)},{1},{},{},options )