% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function phip = barrier_sum( phi1,phi2 )
    phip.fun = @(x)(phi1.fun(x) + phi2.fun(x));
    phip.diff= @(x)(phi1.diff(x) + phi2.diff(x));
    phip.hess= @(x)(phi1.hess(x) + phi2.hess(x));
    phip.nu = phi1.nu + phi2.nu;
    phip.member = @(x)( phi1.member(x) && phi2.member(x));
    phip.activeset = @(x)( phi1.activeset(x) );
end