%% [X1,pval,output] = nlsdp_solve_frank_wolfe(X0,f,Aeq,beq,A,b,options)
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)

function [X1,pval,output] = nlsdp_solve_frank_wolfe(X0,f,Aeq,beq,A,b,options)

    % set default options
    if ~isfield(options,'epsilon')
        options.epsilon = 1e-4;
    end
    if ~isfield(options,'verbose')
        options.verbose = 'iterations';
    end
    
    % initialise printed output
    print_welcome()
    
    % initialise timer
    tic
    
    % find feasible initial point
    d = length(X0);
    cvx_begin sdp quiet
        variable X0(d,d) hermitian
        variable t
        maximise t
        subject to
            X0 >= t*eye(d);
            t >= 0;
            for i = 1:length(A)
                trace(X0*A{i}) <= b(i);
            end
            for i = 1:length(Aeq)
                trace(X0*Aeq{i}) == beq(i);
            end
    cvx_end
    
    % initialise main loop
    X = full(X0);
    pval = f.fun(X);
    eps = options.epsilon;
    eps_iter = Inf;
    iter = 0;
    tstamps = [];
    epsstamps = [];
    
    print_columns( options );
    print_update( iter,pval,[],[],options )
    
    while eps_iter > eps
        [V,dval,y] = FW_iteration(X,f,Aeq,beq,A,b);
        options2 = optimset('TolFun',options.epsilon/2);
        if strcmp( f.conv,'convex' )
            linfun = @(t)(f.fun(X+t*V));
            [t,pval] = fminbnd(linfun,0,1,options2);
        elseif strcmp( f.conv,'concave' )
            linfun = @(t)(-f.fun(X+t*V));
            [t,pval] = fminbnd(linfun,0,1,options2);
            pval = -pval;
        end
        %[t,pval,dval,max(eig(V))]
        pval = real(pval);
        X = X + t*V;
        eps_iter = abs(pval - dval);
        iter = iter + 1;
        print_update( iter,pval,eps_iter,dval,options )
        
        % convergence timing
        t_iter = toc;
        tstamps(end+1) = t_iter;
        epsstamps(end+1) = eps_iter;
    end
    X1 = X;
    output.eps_iter = eps_iter;
    output.lagrangemultiplier = cell2mat(y);
end

function [V,dval,y] = FW_iteration(X,f,Aeq,beq,A,b)
    if strcmp( f.conv,'convex' )
        c = 1;
    elseif strcmp( f.conv,'concave' )
        c = -1;
    end
    X = full((X +X')/2);
    d = length(X);
    grad = f.diff(X);
    grad = (grad + grad')/2;
    cvx_begin sdp quiet
        variable Xp(d,d) hermitian
        dual variable y{length(Aeq)}
        minimize c*trace(Xp*grad)
        subject to
            Xp >= 0;
            for i = 1:length(A)
                trace(Xp*A{i}) <= b(i);
            end
            for i = 1:length(Aeq)
                y{i}: trace(Xp*Aeq{i}) == beq(i);
            end
    cvx_end
    V = value(Xp) - X;
    dval = real(f.fun(X)+c*real(trace(grad*V)));
end
function [] = print_welcome()
    fprintf('\n')
    fprintf('NonlinSDP (C) 2024 Thomas Van Himbeeck\n')
    fprintf('Frank-Wolfe solver\n')
end
function [] = print_columns( options )
    if strcmp(options.verbose,'iterations') 
        fprintf('|iter|   fval   |    eps    |  dual val |\n');
    end
end

function [] = print_update(count,primal,eps,dual,options)
    
    scount = num2str(count,'%4u');
    sprimal = num2str(primal,'%8.5f');
    slambda = num2str(eps,'%7.1e');
    sdual = num2str(dual,'%9.5f');
    if strcmp(options.verbose,'iterations') 
        fprintf('|%4s| %8s |   %7s | %9s |\n',scount,sprimal,slambda,sdual);
    end
end
        