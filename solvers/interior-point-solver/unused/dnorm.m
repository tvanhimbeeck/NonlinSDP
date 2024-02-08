% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function y = dnorm(c,x,F)
    H = F{3}(x);
    y = sqrt(c'*(H\c));
end