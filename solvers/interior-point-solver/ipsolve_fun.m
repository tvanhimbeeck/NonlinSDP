%% [f_opti,x_opti,output] = ipsolve_fun( f,phi,x0,beta,A,b,options )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Computes the minimum 
%       min     f(x) 
%       st.     x \in S
%               Aeq.x = beq
%               Aeq.x < beq
% of the convex function f over a set S, where 
%   - phi is a self-concordant barrier for S with parameter nu
%   - f is beta-compatible with phi
%   - x0 is strictly feasible

function [f_opti,x_opti,output] = ipsolve_fun( f,phi,x0,Aeq,beq,A,b,options )
    
    delta = 0.1; % chosen arbitrarily
    d = length(x0);
    F1 = epigraph_barrier( f,phi,d,f.fun(x0)+2*delta );
    F2 = gen_ineqconsbarrier([A,zeros(size(A,1),1)],b);
    F = barrier_sum(F1,F2);
    z0 = [x0;f.fun(x0)+delta];
    Aeq = [Aeq,zeros(size(Aeq,1),1)];
    if strcmp( f.conv,'convex' )    
        [f_opti,z_opti,output] = ipsolve_set( [zeros(d,1);1],F,z0,Aeq,beq,options);
    elseif strcmp( f.conv,'concave' )
        [f_opti,z_opti,output] = ipsolve_set( [zeros(d,1);-1],F,z0,Aeq,beq,options );
        f_opti = -f_opti;
    end
    x_opti = z_opti(1:d);
end    