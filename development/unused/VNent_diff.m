function dfval = VNent_diff( rho )
    dfval = - logm( rho ) - eye(length(rho));
    dfval = (dfval + dfval')/2/log(2);
end