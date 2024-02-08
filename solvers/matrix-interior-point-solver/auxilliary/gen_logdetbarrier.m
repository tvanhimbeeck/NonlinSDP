% computes the logdet barrier, its gradient and hessian in the Gell-Mann
% basis

function phi = gen_logdetbarrier( D )
    G = gen_GellMann(D);
    phi = struct;
    phi.fun = @(x)( logdet( x,G ) );
    phi.diff = @(x)( logdet_diff( x,G ) );
    phi.hess = @(x)( logdet_hess( x,G ) );
    phi.nu = D;
    phi.member = @(x)( member( x,G ) );
    
end

function y = logdet( x,G )
    y = -log( det( vec2mat( x,G ) ));
end

function g = logdet_diff( x,G )
    g = -mat2vec( inv( vec2mat(x,G) ),G );
end

function H = logdet_hess( x,G )
    l = length( G );
    H = zeros( l );
    drho = decomposition( vec2mat( x,G ),'chol' );
    for i = 1:l
        rhoGrhoi = (drho\G{i})/drho;
        for j = 1:i
            H(j,i) = real(inner_prod( G{j},rhoGrhoi ));
            H(i,j) = H(j,i);
        end
    end
end

function check = member( x,G )
    [~,flag] = chol( vec2mat(x,G) );
    check = ~flag;
    if check == 0
        fprintf('not sdp')
    end
end