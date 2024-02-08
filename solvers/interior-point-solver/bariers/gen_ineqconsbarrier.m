% Copyright (C) 2022 Thomas Van Himbeeck (Licence: GLPv3)
function phi = gen_ineqconsbarrier( A_vec,b_vec )
    phi.fun = @(x)( barrier(x,A_vec,b_vec));
    phi.diff = @(x)( barrier_diff(x,A_vec,b_vec));
    phi.hess = @(x)( barrier_hess(x,A_vec,b_vec));
    phi.nu = size(A_vec,1);
    phi.member = @(x)( member( x,A_vec,b_vec ) );
end

function y = barrier(x,A,b)
    y = 0;
    for i = 1:size(A,1)
        y = y - log( - A(i,:)*x + b(i));
    end
end
function g = barrier_diff(x,A,b)
    g = 0;
    for i = 1 :size(A,1)
        g = g + A(i,:)'/( - A(i,:)*x + b(i));
    end
end
function H = barrier_hess(x,A,b)
    H = 0;
    for i = 1 :size(A,1)
        H = H + A(i,:)'*A(i,:)/( - A(i,:)*x + b(i))^2;
    end
end
function bool = member(x,A,b)
    bool = true;
    for i = 1:size(A,1)
        bool = bool && ( A(i,:)*x<b(i) );
    end
end