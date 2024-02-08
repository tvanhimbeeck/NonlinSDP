%% f2 = compose_with_linfun(f1,A)
% constructs the function f2(rho) := f1(A*rho*A')

function f2 = compose_with_linfun(f1,A)
    f2.fun  = @(rho)    (   f1.fun (A*rho*A')         );
    f2.diff = @(rho)    (A'*f1.diff(A*rho*A')       *A);
    f2.hess = @(rho,X)  (A'*f1.hess(A*rho*A',A*X*A')*A);
    f2.conv = f1.conv;
    f2.beta = f1.beta;
end
