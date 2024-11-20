function y = VNent(X)
    y = - trace(X*logm(X));
end