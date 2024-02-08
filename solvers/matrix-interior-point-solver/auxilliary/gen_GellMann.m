function G = gen_GellMann( D )

    G = cell( 1,D^2 );
    
    G{1} = eye( D )/sqrt(D);
    counter = 2;
    
    for i = 2:D
        mat = diag([ ones(1,i-1),(1-i),zeros(1,D-i) ])/sqrt( (i-1) + (i-1)^2 );
        G{counter} = sparse( mat );
        counter = counter +1;
    end
    for i = 1:D
        for j = i+1:D
            mat = zeros(D);
            mat(i,j) = 1/sqrt(2);
            mat(j,i) = 1/sqrt(2);
            G{counter} = sparse( mat );
            counter = counter +1;
        end
    end
    
    for i = 1:D
        for j = (i+1):D
            mat = zeros(D);
            mat(i,j) = 1i/sqrt(2);
            mat(j,i) = - 1i/sqrt(2);
            G{counter} = sparse( mat );
            counter = counter +1;
        end
    end
end