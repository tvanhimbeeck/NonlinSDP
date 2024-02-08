function H = hess_determinant( rho,rho_det,G )
    % give hessian in basis G{i}
    L = length( G );
    H = zeros( L );
    for i =1:L
        rhoG{i} = rho\G{i};
    end
    for i = 1:L
        for j = 1:L
            H(i,j) = rho_det*real(trace( rhoG{i} ))*real(trace( rhoG{j} ))... 
                        - rho_det*real(trace( (rhoG{i})*(rhoG{j}) ));
        end
    end
    H = real(H + H')/2;
end