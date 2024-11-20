% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function Y = power_grad(p,rho,X)
    % gradient d()^p|_{\rho}[X]
    
    d = length(rho);
    [V,D] = eig(rho,'vector');
    
    [A,B] = meshgrid(D,D);
    Ap = reshape(A,[1,d^2]);
    Bp = reshape(B,[1,d^2]);
    Xp = reshape(V'*X*V,[1,d^2]);
    Xv = V'*X*V;
    
%     factor = sin(pi*p)/pi;
%     integrant = @(s)(s/(1-s))^(p)./(s+Ap*(1-s)).*Xp./(s+Bp*(1-s));
    
    factor2 = sin(pi*p)/(pi*(1-p));
    integrant2 = @(s)s^p*(1-s)^(1-p)*(  1./(s+Ap*(1-s)).*Xp.*Bp./(s+Bp*(1-s)).^2 ...
                                    + Ap./(s+Ap*(1-s)).^2.*Xp./(s+Bp*(1-s))   );
                                
    integrant3 = @(s) (integrant_fun(s,p,D,Xv,d));
                                
    Yp = factor2*integral(integrant3,0,1,'ArrayValued',true,'RelTol',0,'AbsTol',1e-5);
    Y = V*reshape(Yp,[d,d])*V';
end

function y = integrant_fun(s,p,D,Xv,d)
    Y = diag(1./(s+D*(1-s)))*Xv*diag(D./(s+D*(1-s)).^2);
    Y = (Y+Y');
    y = s^p*(1-s)^(1-p)*reshape(Y,[1,d^2]);
end