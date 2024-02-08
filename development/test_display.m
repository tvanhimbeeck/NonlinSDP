

function [] = test_display( options  )
    

    wait = 1;
    % start algorithm
    print_welcome(options)
    
    
    % validate and print initial point
    x0 = rand(4,1);
    c = rand(4,1);
    tau = 0.2;
    count1 = 0;
    lambda = rand(1);
    
    print_initialpoint(c'*x0,lambda,options);

    
    % auxilliary path following scheme
    
    t = 1;
    y = x0;
    
    while lambda > tau
        
        % comput next point
        count1 = count1 + 1;
        ynew = rand(4,1);
        lambda = rand(1);
        pause(wait);
        
        print_auxupdate( count1,c'*ynew,lambda,options )
        
    end
    
        
        
    % intermediate point
        
    count2 = 1;

    xnew = rand(4,1);
    eps = rand(1);
    pause(wait);
    
    print_mainupdate( count1,count2,c'*xnew,eps,c'*xnew-eps,options )
    
    x = xnew;
    
    % main path-following scheme

    
    t=0;
    tmax = 10;
    
    while t < tmax

        % compute new point
        count2 = count2 + 1;
        tnew = t +1;
        xnew = rand(4,1);
        pause(wait);
        x = xnew;
        t = tnew;
        eps_iter = rand(1);
        
        print_mainupdate( count1,count2,c'*x,eps_iter,c'*x-eps_iter,options )

    end
        
    print_ending( count1+count2,options )
end


% functions for printing
function [] = print_welcome(options)

    if strcmp(options.verbose,'iterations')||strcmp(options.verbose,'real-time')
        fprintf('quantum-ip-solver (C) 2023 Thomas Van Himbeeck\n\n')
        print_columns()
        fprintf('\n')
    end
end
function [] = print_initialpoint( fval,lambda,options )
    count1 = 0;
    if strcmp(options.verbose,'iterations')
        status_initialpoint();
        fprintf('\n')
        print_update(count1,fval,lambda,[]);
        
    elseif strcmp(options.verbose,'real-time')
        print_update(count1,fval,lambda,[]);
        fprintf('\n');
        status_initialpoint();
    end
end
function [] = print_auxupdate( count1,fval,lambda,options )
    if strcmp(options.verbose,'iterations')
        if count1==1
            fprintf('\n')
            status_auxpath()
        end
        fprintf('\n')
        print_update(count1,fval,lambda,[]);
    
    elseif strcmp(options.verbose,'real-time')
        erase_line(); % erase status
        fprintf('\b');
        erase_line(); % erase update
        print_update(count1,fval,lambda,[]);
        fprintf('\n');
        status_auxpath();
    end
end
function [] = print_mainupdate( count1,count2,fval,lambda,dval,options )
    if strcmp(options.verbose,'iterations')
        if count2 ==1
            fprintf('\n');
            status_mainpath();
        end
        fprintf('\n');
        print_update(count1+count2,fval,lambda,dval);
    elseif strcmp(options.verbose,'real-time')
        erase_line(); % erase status
        fprintf('\b');
        erase_line(); % erase update
        print_update(count1+count2,fval,lambda,dval);
        fprintf('\n');
        status_mainpath();
    end
end
function [] = print_columns()
    fprintf('|iter|   fval   | eps/lambda|  dual val |')
end
function [] = status_auxpath()
    fprintf('|> Status: finding centering point...   |')
end
function [] = status_mainpath()
    fprintf('|> Status: finding optimum...           |')
end
function [] = status_initialpoint()
    fprintf('|> Status: initial point                |')
end
function [] = erase_line()
    for i=1:41
        fprintf('\b')
    end
end
function [] = print_update( count,primal,lambda,dual )
    
    scount = num2str(count,'%4u');
    sprimal = num2str(primal,'%8.5f');
    slambda = num2str(lambda,'%7.1e');
    sdual = num2str(dual,'%9.5f');
    fprintf('|%4s| %8s |   %7s | %9s |',scount,sprimal,slambda,sdual);
end
function [] = print_ending( count,options )
    if strcmp(options.verbose,'iterations')
        fprintf('\n');
    elseif strcmp(options.verbose,'real-time')
        erase_line(); % erase status
    end
    fprintf('|> Solved: solution found in %4u steps |\n\n',count );
end