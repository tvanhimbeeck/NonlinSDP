% D dimensional QKD
D = 4;
K = 2;
P = 4;
beta = 0.3;
epsilon = 1e-4;

% KrausP{1} = eye( D )*sqrt(p(1));
% KrausP{2} = eye( D )*sqrt(p(2));
% [I,J] = meshgrid(0:(D-1),0:(D-1));
% U = exp(1i*2*pi/D*I.*J)/sqrt(D);
% for k = 1:K
%     KrausKP{k,1} = ((1:D)==k)*sqrt(p(1));
%     KrausKP{k,2} = KrausKP{k,1}*U*sqrt(p(2));
% end

povm = rand_keyMap( K,P,D );

f = gen_keyrate( povm,'povm' );
%f = gen_cond_entropy( povm,'povm' );

X0 = eye( D )/D;
%rho0 = rand_state(D);
Aeq = {}; beq = {};
Aeq{1} = eye(D);
beq{1} = 1;
A = {};
b = {};

%A{1} = - f.diff(rho0);
%b{1} = trace(A{1}*rho0)+1;
% [C,d,rho0] = rand_constraints( 20,D );
% for i = 1:10
%     Aeq{i} = C{i};
%     beq{i} = d{i};
%     A{i} = C{i+10};
%     b{i} = d{i+10}+0.5;
% end

%profile clear
%profile on -timestamp

options.solver = 'frank-wolfe';
tic
[rho1,fval,output] = nlsdp_solve(rho0,f,Aeq,beq,A,b,options);
toc
options.solver = 'interior-point';
tic
[rho1ip,fvalip,outputip] = nlsdp_solve(rho0,f,Aeq,beq,A,b,options);
toc
[fval,fval-output.eps_iter,fvalip,fvalip-outputip.eps_iter]
colortheme = hsv;
semilogy(tstamps,epsstamps,'Color',colortheme(floor(263*count/5)+1,:))
hold on
semilogy(tstampsip,epsstampsip,'Color',colortheme(floor(263*count/5)+1,:),'Linestyle','--')
count = count+1;
%profile viewer
%profile off
