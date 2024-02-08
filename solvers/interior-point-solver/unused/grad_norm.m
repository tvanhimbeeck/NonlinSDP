% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function y = grad_norm(x,F)
    y = dnorm(F{2}(x),x,F);
end