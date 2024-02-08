function M = vec2mat( x,G )
    D = sqrt(length(x));
    M = full(zeros( D ));
    for i = 1:D^2
        M = M + x(i)*G{i};
    end
    M = (M+M')/2;
end