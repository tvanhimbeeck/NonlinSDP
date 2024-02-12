%% q = eval_povm( rho,povm )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)

function q = eval_povm( rho,povm ) 

    for i = 1:length(povm)
        q(i) = trace( povm{i}*rho );
    end
end