% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function Y = power_hess(p,rho,X,A)

    % hessian d^2()^p|_rho[X,A]
    
    d = length(rho);
    [V,D] = eig(rho,'vector');
    
    Xp = V'*X*V;
    Ap = V'*A*V;
    
    factor = sin(pi*p)/(pi);
    integrant = @(s) (integrant_fun(s,p,D,Ap,Xp,d));
                                
    Yp = factor*integral(integrant,0,1,'ArrayValued',true,'RelTol',0,'AbsTol',1e-5);
    Y = V*reshape(Yp,[d,d])*V';
end

function y = integrant_fun(s,p,D,Ap,Xp,d)
    Y = diag(1./(s+D*(1-s)))*Xp*diag(1./(s+D*(1-s)))*Ap*diag(1./(s+D*(1-s)));
    Y = (Y+Y');
    y = -s^p*(1-s)^(1-p)*reshape(Y,[1,d^2]);
end