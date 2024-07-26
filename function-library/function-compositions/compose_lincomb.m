%% f2 = compose_lincomb( f,a,K )
% Constructs the matrix function 
%   	f2(X) := sum_i a(i) * f{i}( K{i}*X*K{i}' ) 
% if the convexities match.
% 
% Variantions with 1 or 2 inputs
%   compose_lincomb( f ) assumes a(i) = 1, K{i} = 1
%   compose_lincomb( f,a ) assumes K{i} = 1

function f2 = compose_lincomb( varargin )
    f = varargin{1};
    if nargin >= 2
        a = varargin{2};
    else
        a = repmat([1],1,length(f));
    end
    if nargin == 3
        K = varargin{3};
    else
        K = repmat({1},1,length(f));
    end
    
    f2.fun  = @(X)    ( fun(f,a,K,X) );
    f2.diff = @(X)    ( diff(f,a,K,X) );
    f2.hess = @(X,V)  ( hess(f,a,K,X,V) );
    f2.beta = max_beta( f );
    f2.conv = convexity_sum( f,a );
    
end

function fval = fun( f,a,K,X )
    fval = 0;
    for i = 1:length(f)
        fval = fval + a(i)*f{i}.fun( K{i}*X*K{i}' );
    end
end

function grad = diff( f,a,K,X )
    grad = 0;
    for i = 1:length(f)
        grad = grad + a(i)*K{i}'*f{i}.diff( K{i}*X*K{i}' )*K{i};
    end
end

function HX = hess( f,a,K,X,V )
    HX = 0;
    for i = 1:length(f)
        HX = HX + a(i)*K{i}'*f{i}.hess( K{i}*X*K{i}',K{i}*V*K{i}' )*K{i};
    end
end

function beta = max_beta( f )
    for i = 1:length(f)
        betas(i) = f{i}.beta;
    end
    beta = max(betas);
end

function conv = convexity_sum( f,a )
    
    % check matching convexities
    for i = 1:length(f) 
        if strcmp(f{i}.conv,'concave')
            c(i) = -sign(a(i));
        elseif strcmp(f{i}.conv,'convex')
            c(i) = sign(a(i));
        end
    end
    if c >= ones(1,length(f))
        conv = 'convex';
    elseif c == -ones(1,length(f))
        conv = 'concave';
    else
        fprintf( 'forbidden composition\n' )
    end
end
    