% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function [A,b,rho] = rand_constraints( n_cons,D )

    A = cell(1,n_cons);
    b = cell(1,n_cons);
    A{1} = eye( D );
    b{1} = 1;
    
    rho = rand_state( D );
    
    for i = 2:n_cons
        mat = rand( D );
        A{i} = mat+mat';
        b{i} = trace( rho*A{i} );
    end
end
