function f_vec = vectorize_function( f,G )
    D = length(G{1});
    f_vec.fun  = @(z)( f.fun( vec2mat( z(1:D^2),G ) ) );
    f_vec.diff = @(z)( mat2vec( f.diff( vec2mat( z(1:D^2),G ) ),G ) );
    f_vec.hess = @(z)( vectorise_map( @(V)(f.hess(vec2mat(z(1:D^2),G),V)),G ) );
    f_vec.conv = f.conv;
    f_vec.beta = max(f.beta);
end

function H = vectorise_map( map,G )
    d = length(G);
    for i = 1:d
        HGi = map(G{i});
        for j = 1:d
            H(i,j) = real(inner_prod(G{j},HGi));
        end
    end
    H = (H + H')/2;
    H = real(H);
end