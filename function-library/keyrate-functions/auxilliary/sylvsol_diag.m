function X = sylvsol_diag(A,B,C)
%SYLVSOL_DIAG  Solve Sylvester equation.
%   X = SYLVSOL_DIAG(A, B, C) is the solution to the Sylvester equation
%   AX + XB = C, where A and B are diagonal

[m,m] = size(A);
[n,n] = size(B);
X = C ./(diag(A)*ones(1,n) + ones(m,1)*diag(B)');

end
