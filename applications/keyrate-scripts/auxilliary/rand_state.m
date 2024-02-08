% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function M = rand_pos( D )
    mat = rand(D);
    M = mat*mat'/trace(mat*mat');
end