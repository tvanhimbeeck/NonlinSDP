%% [rho1,fval,output] = ipsolve_matrixfun( X0,fun,Aeq,beq,A,b,options )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)

function [X1,fval,output] = ipsolve_matrixfun( X0,f,Aeq,beq,A,b,options )
    
    % Default options
        if ~isfield(options,'domain')
            options.domain = 'complex';
        end
    
    % Protocol parameters
        D = length( X0 ); % dimension of the problem
        n_eqcons = length( Aeq ); % number of equality constraints
        n_cons = length( A ); % number of inequality constraints
    
    % Determine basis
        if strcmp(options.domain,'complex')
            G = gen_GellMann( D );
        end
        
    % Vectorize everything in Gell-Mann basis

        x0 = mat2vec( X0,G ); 
        
        % equality constraints
        Aeq_vec = [];
        beq_vec = [];
        for i = 1:n_eqcons
            Aeq_vec(i,:) = mat2vec( Aeq{i},G )';
            beq_vec(i,1) = beq(i);
        end
        
        % inequality constraints
        A_vec = [];
        b_vec = [];
        for i = 1:n_cons
            A_vec(i,:) = mat2vec( A{i},G )';
            b_vec(i,1) = b(i);
        end
        
        % can be replaced with f_vec = vectorize_function( f,G );
        
        % objective functions
        f_vec.fun  = @(x)( f.fun( vec2mat( x,G ) ) );
        f_vec.diff = @(x)( mat2vec( f.diff( vec2mat( x,G ) ),G ) );
        f_vec.hess = @(x)( vectorise_map( @(V)(f.hess(vec2mat(x,G),V)),G ) );
        f_vec.conv = f.conv;
        f_vec.beta = f.beta;
        
    % Optimisation
        phi = gen_logdetbarrier( D );
        
        [fval,x1,output] = ipsolve_fun( f_vec,phi,x0,Aeq_vec,beq_vec,A_vec,b_vec,options );
        X1 = vec2mat( x1,G );

end

function H = vectorise_map( map,G )
    d = length(G);
    for i = 1:d
        HGi = map(G{i});
        for j = 1:d
            H(i,j) = real(inner_prod(G{j},HGi));
        end
    end
    H = (H + H')/2;
    H = real(H);
end

        