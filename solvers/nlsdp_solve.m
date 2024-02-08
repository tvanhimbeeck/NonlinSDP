%% [X1,fval,output] = nlsdp_solve( X0,fun,Aeq,beq,A,b,options )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)

% Solves convex nonlinear sdp problems optimisation problems of the
% form:
%           minimize      fun( X )
%           subject to    Tr[ Aeq{i} X ] == beq{i}
%                         Tr[ A{i}   X ] <= b{i}
%                         X >= 0
% where fun(X) is a convex differentiable function defined for hermitian 
% positive semidefinite operator (or the corresponding maximisation problem for
% concave functions).
%
% Parameters:
%       - X0  = starting point (d-by-d strictly feasible hermitian matrix)
%       - fun = objective function satisfying the function structure template
%       - Aeq, A = cell arrays with d-by-d hemitian matrices
%       - beq, b = matching cell arrays of numbers
%
%       - options.algorithm = 'frank-wolfe' or 'interior-point'

function [X1,fval,output] = nlsdp_solve( X0,fun,Aeq,beq,A,b,options )
    
    if ~isfield(options,'solver') || strcmp( options.solver,'interior-point' )
        [X1,fval,output] = ipsolve_matrixfun( X0,fun,Aeq,beq,A,b,options );
    
    elseif strcmp( options.solver,'frank-wolfe' )
        [X1,fval,output] = nlsdp_solve_frank_wolfe( X0,fun,Aeq,beq,A,b,options );      
    
    end
    
end