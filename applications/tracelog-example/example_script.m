d = 2;

X0 = rand_state(2);
Aeq = {eye(d),[0 1i;-1i 0]};
beq = {1,0};
A = { rand_herm(2),rand_herm(2),rand_herm(2)};
b = { trace(X0*A{1})+0.01,trace(X0*A{2})+0.01,trace(X0*A{3})+0.01 };

options.solver = 'interior-point';
nlsdp_solve( X0,gen_trlog,Aeq,beq,A,b,options )

function X = rand_state(d)
    X = rand(d);
    X = X*X'/trace(X*X');
end
function X = rand_herm(d)
    X = rand(d)-0.5*ones(d);
    X = X+X';
end