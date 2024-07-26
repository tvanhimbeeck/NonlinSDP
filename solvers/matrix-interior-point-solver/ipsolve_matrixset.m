%% [rho1,fval,output] = ipsolve_matrixset(X0,f,g,F,Aeq,beq,options )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)
%
% Solves the problem
%       minimize    f(X)+g(y)
%       such that   F(X)<y
%                   X>= 0
%                   tr(Aeq{i}*X) = beq{i}
%
% where f,F are convex matrix function and g is a convex scalar function

function [X1,fval,output] = ipsolve_matrixset( X0,f,g,F,Aeq,beq,options )
    
    % Protocol parameters
        delta = 0.01;
        D = length( X0 ); % dimension of the problem
        n_eqcons = length( Aeq ); % number of equality constraints
    
    % Vectorize everything in Gell-Mann basis
        G = gen_GellMann(D);
        x0 = mat2vec( X0,G );
        y0 = F.fun( X0 ) + delta;
        z0 = [x0;y0];
        
        % equality constraints
        Aeq_vec = [];
        beq_vec = [];
        for i = 1:n_eqcons
            Aeq_vec(i,:) = [mat2vec( Aeq{i},G )',0];
            beq_vec(i,1) = beq(i);
        end
        

        % objective function
        f_vec.fun  = @(z)( f.fun( vec2mat( z(1:D^2),G ) ) + g.fun( z(D^2+1) ) );
        f_vec.diff = @(z)( [mat2vec( f.diff( vec2mat( z(1:D^2),G ) ),G );...
                            g.diff( z(D^2+1) ) ] );
        f_vec.hess = @(z)( [vectorise_map( @(V)(f.hess(vec2mat(z(1:D^2),G),V)),G ),zeros(D^2,1);...
                            zeros(1,D^2), g.hess(z(D^2+1)) ] );
        f_vec.conv = f.conv;
        f_vec.beta = max(f.beta,g.beta);
        
        % epigraph
        phi = gen_logdetbarrier( D ); % SDP cone of dimension d^2
        F_vec = vectorize_function( F,G );
        phi2 = epigraph_barrier( F_vec,phi,D^2,0 ); % epigraph of dimension d^2+1 
        
    % Optimisation
        
        [fval,z1,output] = ipsolve_fun( f_vec,phi2,z0,Aeq_vec,beq_vec,{},[],options );
        X1 = vec2mat( z1(1:D^2),G );

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
        