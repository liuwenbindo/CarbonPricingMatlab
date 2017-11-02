function [ f , x , fcount , gcount,iter ] = quasi_newton( myfun , x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           CAYUGA RESEARCH, Oct 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quasi-Newton method (Algorithm 6.1) described in the book of chapter 6.
%
% Input:
% myfun - function's name
% x0 - start point
% DELETED varargin - additional parameters for myfun
%
% Output:
% f - final function value
% x - final point where we stop
% fcount - number of evaluation of function value
% gcount - number of evaluation for first gradient of f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol=1e-4; % changed from 1e-4 to 1e-2 due to the computing time
itbnd=500*1e20;
n=length(x0);
H=eye(n);
x=x0;
[f,g]=feval(myfun,x);% x: ROW, g: ROW
g = g'; % g: COL
fcount=1;
gcount=1;
iter=0;
grads = 0;

while(norm(g)>tol && iter<itbnd && fcount<1e3 ) % changed fcount < 1e3 to 1e2 due to the computing time 
    p=-H*g; % H: n*n, g: n*1 ==> p: n*1
    if (g'*p>0) % (1*n) * (n*1)
        H=eye(n);
        p=-H*g;
    end
    [alpha,fc,gc]=line_search(myfun,f,g,x,p');
    x=x+alpha*p';
    s=alpha*p;
    [f,g_new]=feval(myfun,x);
    
    fcount=fcount+fc+1;
    gcount=gcount+gc+1;
    y=g_new-g';
    y = y';
    g=g_new;    
    g = g';
    
    if (iter==0)
        H=y'*s/(y'*y)*eye(n);
    end
     
    iter=iter+1;
    rho=1/(y'*s);
    %H=H-(H*y*y'*H)/(y'*H*y)+(s*s')/(y'*s);
    H=(eye(n)-rho*s*y')*H*(eye(n)-rho*y*s')+rho*s*s';
    fprintf('Now we are running No. ');
    iter
    fprintf(' th Loop. Now the norm of gradient is ');
    norm(g)
    fprintf('\n');
    grads = [grads,norm(g)];
end

end