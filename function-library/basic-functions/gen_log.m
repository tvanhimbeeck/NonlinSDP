%% f = gen_log()

function f = gen_log()
        
    f.fun = @(x)   ( log(x) );
    f.diff = @(x)  ( 1/x );
    f.hess = @(x)  ( -1/x^2 );
    f.conv = 'concave';
    f.beta = 2;
    f.input = 'scalar';
end