%% povm = rand_povm( K,P,D )
% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
%
% Generates a random POVM indexed by two indices summing to identity

function povm = rand_povm( K,P,D )
    
    emptyCell = cell(K,P); 
    povm = cellfun( @(x)rand_pos( D ),emptyCell,'UniformOutput',false );
    
    Sum = 0;
    for i = 1:numel( povm )
        Sum = Sum + povm{i};
    end
    K = inv(sqrtm(Sum));
    povm = cellfun( @(x)K*x*K,povm,'UniformOutput',false );
end