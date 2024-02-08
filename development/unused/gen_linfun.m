function gen_linfun( A )
    % implements the linear matrix function Tr[A*X]
G = gen_GellMann(length(A));
linfun = struct;
linfun.fun  = @(x) ( trace( A*vec2mat(x,G) ) );
linfun.diff = @(x) ( mat2vec( A,G ) );
linfun.hess = @(x) ( zeros(length(A)) );
linfun.beta = 1;
end