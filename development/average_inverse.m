function lambda = average_inverse(H,v)
    K = chol(H);
    y = K'\v;
    lambda = y'*y;
end