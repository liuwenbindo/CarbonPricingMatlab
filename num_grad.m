% [df] = num_grad(func,X)
%
% Function to compute the numerical gradient of an arbitrary objective
% function.
%
% Inputs:
% func: Function handle for which numerical derivative is to 
%		be obtained.
% X:    Point of interest about which derivative is to be
%		obtained.
% Outputs:
% df:   Numerical derivative of function func (vector of size
%		n=length(X)).


function df = num_grad(func, X)

df = zeros(length(X), 1);
h = 1e-5;
% for each dimension of objective function
for i=1:length(X)
    x1 = X;
    x2 = X;
    % vary variable i by a small amount (left and right)
    x1(i) = X(i) - h;
    x2(i) = X(i) + h;
    % evaluate the objective function at the left and right points
    y1 = func(x1);
    y2 = func(x2);
    % calculate the slope (rise/run) for dimension i  
    df(i) = (y2 - y1) / (2*h);
end

end
