m=100;
myfun = @utilitycalc;
dim = 63;
x = double(zeros(m,63));
for ii =1:m
    x(ii,:) = 2 * rand(1,63);
end
gamma=1; % a parameter in RBF (phi function)
deg=-1; % in the RBF model, p(x) is set to zero.
% evaluate the objective fcn at each of these m points
f=zeros(m,1);
for jj=1:m
    xjj_size = size(x(jj,:));
    input = x(jj,:);
    f(jj)= utilitycalc(input);
end

myfun = ADfun('utilitycalc',1);
ind=find(f==min(f));
%
% set the mono decreasing lambda sequence
Ls=[10*0.3.^(0:2),0];
%
% set the global minimizer guesstimate
xstar=2*(2*rand(dim,1)-1);
%
% gradient used?
useg=1; % gradient is not used
%
% set the function used in the RBF approximation
method='cubic';

[fmin1,xmin1,iter1,xbar_min1,fbar_min1,fcount]=global_RBF_TRM(myfun,x,f,ind,xstar,Ls,method,deg,gamma,useg);


