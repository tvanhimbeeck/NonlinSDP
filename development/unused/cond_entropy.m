function fval = cond_entropy( rho,keyMap )

    [K,P] = size( keyMap );
    [KrausKP,KrausP] = gen_kraus( keyMap );
    
    % % Perturb matrix to avoid structly negative eigenvalues
    % 
    % eigMin = lambda_min(rho);  % check the minimum eigenvalue of this density matrix
    % dim = size(rho,1);
    % if eigMin <= 0
    %     epsilon = -eigMin; % choose the epsilon value for perturbation
    %     rho = (1-epsilon)*rho + epsilon*eye(dim)/dim;
    % end

    % Function evaluation

    fval = 0;
    for p = 1:P
        fval = fval - VNent( KrausP{p}*rho*KrausP{p} );
        for k = 1:K
            fval = fval + VNent( KrausKP{k,p}*rho*KrausKP{k,p} );
        end
    end
    fval = real(fval);

end