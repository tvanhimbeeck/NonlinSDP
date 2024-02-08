% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function F = epigraph_barrier( f,phi,d,u )
    F.fun = @(x)( -log( -f.fun(x(1:d)) + x(d+1) )... 
                 + f.beta^3*phi.fun(x(1:d))... 
                 - log(u-x(d+1))...
                );
    
    F.diff = @(x)(  [ f.diff(x(1:d))           ;-1       ]/(-f.fun(x(1:d)) + x(d+1))...
                 + [ f.beta^3*phi.diff(x(1:d))   ;0        ]...
                 + [ zeros(d,1)         ;1/(u-x(d+1))]...
               );
    F.hess = @(x)( [ f.hess(x(1:d)),    zeros(d,1);
                     zeros(1,d), 0          ]/(-f.fun(x(1:d)) + x(d+1)) ...
                 + [ -f.diff(x(1:d));1]*[-f.diff(x(1:d))',1]./(-f.fun(x(1:d)) + x(d+1))^2 ...
                 + [ f.beta^3*phi.hess(x(1:d)),  zeros(d,1);
                     zeros(1,d), 1/(u-x(d+1))^2]...%%mod
               );
    F.nu = 1 + f.beta^3 *phi.nu + 1;
    F.member = @(x)(memberfun(x,f,phi,d,u));
    F.activeset = @(x) [f.diff(x(1:d));-1]/(-f.fun(x(1:d)) + x(d+1));
end

function bool = memberfun(x,f,phi,d,u)
    if f.fun(x(1:d)) >= x(d+1)
        %fprintf('epigraph problem')
        %disp(f.fun(x(1:d))- x(d+1))
    end
    bool = (f.fun(x(1:d)) < x(d+1)) && phi.member(x(1:d)) && (x(d+1)<u);
end