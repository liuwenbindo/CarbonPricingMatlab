% sacle-valued function tildeFw:
%
%        tildeFw = [y-tildeF(x(1:n))]' * w
%
function fval = tildeFw(x, Extra)

% get the initial values of function.
n = length(x);
y = Extra.y;
w = Extra.w;

fvec = zeros(n,1);
i = 2:(n-1);
fvec(i) = (3-2*x(i)).*x(i)-x(i-1)-2*x(i+1)+ones(n-2,1);
fvec(n) = (3-2*x(n)).*x(n)-x(n-1)+1;
fvec(1) = (3-2*x(1)).*x(1)-2*x(2)+1;

fvec = y- fvec;

fval = fvec'*w;


