function f = gen_VNent() 

    load tau_r8_zeros
    
    f.fun = @(rho)   ( VNent_fun(      rho ) );
    f.diff = @(rho)  ( VNent_diff( rho ) );
    f.hess = @(rho,X)( VNent_hess( rho,pzero,qzero ) );
    f.conv = 'concave';
    f.beta = 1;
end

function fval = VNent_fun(rho)
    e = eig(rho);
    fval = -e'*log2(e+(e==0));
end

function grad = VNent_diff(rho)
    grad = - logm(X)- eye(lenght(rho));
end

function HX = VNent(rho,X)
    HX = -logm_frechet_pade_herm( rho,X,pzero,qzero );
end