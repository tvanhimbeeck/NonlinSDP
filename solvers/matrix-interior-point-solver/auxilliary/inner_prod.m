%% Inner product Tr[A'*B] between two matrices 
% optimized for when A is sparse
function y = inner_prod(A,B)
    if issparse(A)
        [i,j,v] = find(A);
        y = 0;
        for k = 1:length(v)
            y = y  + v(k)'*B(i(k),j(k));
        end
    else
        y = trace(A*B);
    end
end