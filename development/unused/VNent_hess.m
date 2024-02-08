function H = VNent_hess( rho,G,K )
    % computes the hessian of S( Kraus*rho*Kraus' ) in rho
    % in the basis G
    
    dim = length(rho);
    rho = (1 - 1e-10)*rho + 1e-10*eye(length(rho));
    state = K*rho*K';
    inv_state = inv(state);
    
    M = inv_o_adj( G,rho,inv_state );
    M2 = left_action( G,inv_state );
    M3 = left_right_action( G,K );
    
    [U,D] = eig(M,'vector');
    X = U*diag(h(D))*U';
    %X = -M\logm( eye(length(M)) - M );
    %X = (X+X')/2;
    H = - M3'*M2*( X )*M3;
    
    H = real(H+H')/2/log(2);
   
end

function y = h(x)
    for i = 1:length(x)
        if abs(x(i))<1e-3
            y(i) = 1 + x(i)/2 + x(i)^2/3;
        else
            y(i) = -log(1-x(i))./x(i);
        end
    end
end

function M = inv_o_adj( G,A,invA )
    L = length( G );
    % assume G is orthonormal basis for Hilbert Schmidt inner product
    % G{i} is hermitian matrix of same dimension as M
    
    M = zeros( L );
    for i = 1:L
        for j = 1:L
            M(i,j) = (i==j) - trace( G{i}*invA*G{j}*A );
        end
    end
    M = (M + M')/2;
end

function M = left_action( G,invA )
    L = length( G );
    M = zeros( L );
    for i = 1:L
        for j = 1:L
            M(i,j) = trace( G{i}*invA*G{j} ) ;
        end
    end
    M = (M + M')/2;
end

function M = left_right_action( G,K )
    L = length( G );
    M = zeros( L );
    for i = 1:L
        for j = 1:L
            M(i,j) = trace( G{i}*K*G{j}*K' ) ;
        end
    end
end

