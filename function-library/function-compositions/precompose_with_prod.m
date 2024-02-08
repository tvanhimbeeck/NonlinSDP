%% f2 = precompose_with_prod(f1,a)
% constructs the function f2(a) := f1(a*x)

function f2 = precompose_with_prod(f1,a)
    f2.fun  = @(x)    (    f1.fun (a*x) );
    f2.diff = @(x)  (  a*f1.diff(a*x) );
    f2.hess = @(x)(a^2*f1.hess(a*x) );
    f2.conv = f1.conv;
    f2.beta = f1.beta; % maybe this changes with the norm of a actually
end
