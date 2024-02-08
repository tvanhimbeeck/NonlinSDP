%% [f_opt,x_opt,output] = ipsolve_set( c,F,x0,A,b,options )
% Copyright (C) 2023 Thomas Van Himbeeck (Licence: GLPv3)
%
% Computes the minimum 
%       min     <c,x> 
%       st.     x \in S
%               A*x = b
% of the linear function <c,x> over a set S, where 
%   - F is a self-concordant barrier for S with parameter nu
%   - A and b parametrise linear equality constraints
%   - x0 is strictly feasible

function [obj,x1,output] = ipsolve_set( c,F,x0,A,b,options )
        
    % Solver default options
    
    options = set_defaults(options);
    
    
    % Convergence timing
    tic
    output.tstamps = [];
    output.epsstamps = [];
    
    % Disable warnings 
    if strcmp(options.warning,'on')
        warning('on','MATLAB:nearlySingularMatrix')
    elseif strcmp(options.warning,'off')
        warning('off','MATLAB:nearlySingularMatrix')
    end
    
    % start algorithm
    print_welcome( options );
    
    % validate and print initial point
    count1 = 0;
    g = F.diff(x0);
    H = F.hess(x0);
    v = F.activeset(x0);    
    db = A*x0 - b;
    [~,lambda] = newton_stepC( g,H,A,db,v,options );
    if checks( F,x0,A,b,lambda,count1,options,'Initial point' ) < 0
        x1 = x0; obj = c'*x0; options.flag = -1; return;
    end
    
    print_initialpoint( c'*x0,lambda,options );

        
    % auxilliary path following scheme

    
    t = 1;
    y = x0;
    df0 = F.diff(x0);
    
    while lambda > options.tau
        
        % comput next point
        count1 = count1 + 1;
        [t,ynew] = path_follow_step( t,y,-df0,g,H,-options.gamma,A,db,v,options );
        
        % validate and analyse point
        db = A*ynew - b;
        g = F.diff(ynew);
        H = F.hess(ynew);
        
        [~,lambda] = newton_stepC( g,H,A,db,v,options );
        [~,beta_act] = newton_stepC( -t*df0 + g,H,A,db,v,options );
        
        if checks( F,ynew,A,b,beta_act,count1,options,'Iteration' )<0
            x1 = y; obj = c'*y; options.flag = -1; return;
        end
        
        print_auxupdate( count1,c'*ynew,lambda,options );
        y = ynew;
    end
    
        
        
    % intermediate point
        
    count2 = 1;
    x = y;

    xnew = y - newton_stepC( g,H,A,db,v,options );
    
    if checks( F,xnew,A,b,0,count1+count2,options,'Initial point' )<0
        x1 = x; obj = c'*x; options.flag = -1; return;
    end
    
    print_mainupdate( count1,count2,c'*xnew,eps,c'*xnew-eps,options );
    x = xnew;
    
    % main path-following scheme
    db = A*x - b;
    g = F.diff(x);
    H = F.hess(x);
    v = F.activeset(x);
    
    t=0;
    tmax = 1/options.epsilon*(F.nu + (options.beta + sqrt(F.nu))*options.beta/(1-options.beta));
    
    while t < tmax

        % compute new point
        count2 = count2 + 1;
        
        [tnew,xnew] = path_follow_step( t,x,c,g,H,options.gamma,A,db,v,options );
        
        bflag = 1;
        count = 1;
        while bflag ==1
            
            bflag = basic_checks( xnew,F,'Iteration' );
            if bflag == -1 || count>options.reiteratemax
                obj = c'*x; x1 = x; options.flag = -1; return;
            elseif bflag == -2
                xnew = (xnew-x)/100 + x;
            else
                break;
            end
            count = count+1;
        end
        
        % analyse and validate new point 
        g = F.diff(xnew);
        H = F.hess(xnew);
        db = A*xnew - b;
        if max(eig(H)<0)==1
            1;
        end
        
        [~,beta_act] = newton_stepC( -t*df0 + g,H,A,db,v,options );
        check = checks( F,xnew,A,b,beta_act,count1+count2,options,'Iteration' );
        if check == -1
            obj = c'*x; x1 = x; options.flag = -1; return;
        end
        
        x = xnew;
        t = tnew;
        eps_iter = 1/t*(F.nu + (options.beta + sqrt(F.nu))*options.beta/(1-options.beta));
                
        print_mainupdate( count1,count2,c'*x,eps_iter,c'*x-eps_iter,options );
        
        output.epsstamps(end+1) = eps_iter;
        output.tstamps(end+1) = toc;

    end

    print_ending( count1+count2,options );
    
    % return values
    obj = c'*x;
    x1 = x;
    output.eps_iter = eps_iter;
    output.flag = 1;
    output.tstamps;
    output.epsstamps;
end

function [t2,x2,lambda_act] = path_follow_step( t1,x1,c,g,H,gamma_step,A,db,v,options )
    
    [~,lambdac] = newton_stepC( c,H,A,db,v,options );
    t2 = t1 + gamma_step/lambdac;
    [delta,lambda_act] = newton_stepC( t2*c + g,H,A,db,v,options );
    x2 = x1 + delta;
end

function flag = checks( F,x,A,b,beta_act,count,options,type )
    
    if ~isreal(x) || max(~isfinite(x))
        fprintf('\n%s is invalid\n',type);
        flag = -1; return;
        
    elseif ~F.member(x)
        fprintf('\n%s is not strictly feasible\n',type);
        flag = -2; return;
        
    elseif  max(abs(A*x - b))>1e-9
        fprintf('\n%s violates equality constraints\n',type);
        flag = -1; return;
        
    elseif beta_act>options.beta
        fprintf('\n%s left central path\n');
        flag = -1;
        
    elseif count > options.maxiter
        fprintf('\nMax iter reached\n')
        flag = -1; return;
    end
    flag = 1;
end
function flag = basic_checks( x,F,type )
    
    if ~isreal(x) || max(~isfinite(x))
        fprintf('\n%s is invalid\n',type);
        flag = -1; return;
        
    elseif ~F.member(x)
        %fprintf('\n%s is not strictly feasible\n',type);
        flag = -2; return;
    end
    flag = 1;
end
function [delta,lambda]= newton_stepC( g,H,A,db,v,options )
    
    d = length(g);
    n = size(A,1);
    
    if strcmp(options.precondition,'on')
        P = eye(d) -(1- 1/(1+sqrt(v'*v)))*(v*v')/(v'*v);
        M = [P*H*P,P*A';A*P,zeros(n)];
        z = M\[-P*g;-db];
        delta = P*z(1:d);
    else
        M = [H,A';A,zeros(n)];
        z = Mold\[-g;-db];
        delta = zold(1:d);
    end
    mu = z(d+1:d+n);
    lambda = sqrt(z'*M*z);
    fun = @(x)(1/(1+(x^2/(1+x))));
    delta = delta * fun(lambda);
end
function options = set_defaults( options )
    if  ~isfield(options,'epsilon')
        options.epsilon = 1e-4;
    end
    if  ~isfield(options,'maxiter')
        options.maxiter = 1e3;
    end
    if ~isfield(options,'verbose')
        options.verbose = 'real-time';
    end
    if ~isfield(options,'warning')
        options.warning = 'off';
    end
    if ~isfield(options,'reiteratemax')
        options.reiteratemax = 2;
    end
    if ~isfield(options,'precondition')
        options.precondition = 'on';
    end
        
    options.constraint_tol = 1e-6;
    options.beta = 0.1262;
    options.tau = 0.29;
    options.gamma = (options.tau - options.beta);
    
end
% functions for printing
function [] = print_welcome(options)

    if strcmp(options.verbose,'iterations')||strcmp(options.verbose,'real-time')
        fprintf('slsdp-solver (C) 2023 Thomas Van Himbeeck');
        fprintf('\n\n');
        fprintf('|iter|   fval   | eps/lambda|  dual val |');
        %fprintf('\n')
    end
end
function [] = print_initialpoint( fval,lambda,options )
    count1 = 0;
    status = '|> Status: initial point                |';
    if strcmp(options.verbose,'iterations')
        fprintf(status);
        fprintf('\n');
        print_update(count1,fval,lambda,[]);
        
    elseif strcmp(options.verbose,'real-time')
        print_update(count1,fval,lambda,[]);
        fprintf('\n');
        fprintf(status);
    end
end
function [] = print_auxupdate( count1,fval,lambda,options )
    status = '|> Status: finding centering point...   |';
    if strcmp(options.verbose,'iterations')
        if count1==1
            fprintf('\n');
            fprintf(status);
        end
        fprintf('\n');
        print_update(count1,fval,lambda,[]);
    
    elseif strcmp(options.verbose,'real-time')
        erase_line(); % erase status
        fprintf('\b');
        erase_line(); % erase update
        print_update(count1,fval,lambda,[]);
        fprintf('\n');
        fprintf(status);
    end
end
function [] = print_mainupdate( count1,count2,fval,lambda,dval,options )
    status = '|> Status: finding optimum...           |';
    if strcmp(options.verbose,'iterations')
        if count2 ==1
            fprintf('\n');
            fprintf(status);
            fprintf('\n');
        end
        %fprintf('\n');
        print_update(count1+count2,fval,lambda,dval);
    elseif strcmp(options.verbose,'real-time')
        if count2 ==1
            fprintf('\n');
        end
        erase_line(); % erase status
        fprintf('\b\b');
        erase_line(); % erase update
        print_update(count1+count2,fval,lambda,dval);
        fprintf('\n');
        fprintf(status);
    end
end
function [] = print_ending( count,options )
    if strcmp(options.verbose,'iterations')
        fprintf('\n');
    elseif strcmp(options.verbose,'real-time')
        erase_line(); % erase status
    end
    fprintf('|> Solved: solution found in %4u steps |\n\n',count );
end

function [] = print_update( count,primal,lambda,dual )
    
    scount = num2str(count,'%4u');
    sprimal = num2str(primal,'%8.5f');
    slambda = num2str(lambda,'%7.1e');
    sdual = num2str(dual,'%8.5f');
    fprintf('|%4s| %8s |   %7s | %9s |',scount,sprimal,slambda,sdual);
end
function [] = erase_line()
    for i=1:41
        fprintf('\b')
    end
end

% function [] = print_welcome()
%     fprintf('quantum-ip-solver (C) 2023 Thomas Van Himbeeck\n')
% end
% function [] = print_columns()
%     fprintf('\n|iter|   fval   | eps/lambda|  dual val |\n')
% end
% function [] = print_beginauxpath()
%     fprintf('\n|Finding centering point...             |\n')
% end
% function [] = print_beginmainpath()
%     fprintf('\n|Finding optimum...                     |\n')
% end
% function [] = erase_update()
%     for i=1:34
%         fprintf('\b')
%     end
% end
% function [] = print_update( count,primal,lambda,dual )
%     
%     scount = num2str(count,'%4u');
%     sprimal = num2str(primal,'%8.5f');
%     slambda = num2str(lambda,'%7.1e');
%     sdual = num2str(dual,'%9.5f');
%     fprintf('|%4s| %8s |   %7s | %9s |\n',scount,sprimal,slambda,sdual);
% end
% function [] = print_ending( count )
%     fprintf('\nSolution found in %4u steps.\n',count);
% end