function x = mat2vec( M,G )
    D = length( M );
    x = zeros(D^2,1);
    for i = 1:D^2
        x(i) = inner_prod(G{i},M);
    end
    x = real(x);
end