myfun = ADfun('utilitycalc',1);
dim = 63;

m = [0.833862047445461,0.986771045676486,0.844791260010868,1.29452599908852,1.02337420555469,1.03977841748325,0.909265854132061,1.10715414356476,1.07129074654436,1.17879428820710,1.13498630255108,1.20495671381198,1.16330210508919,1.21134083349312,0.916642214563293,1.02588970956374,0.989552264507359,1.02911933651639,0.989687505001852,1.03469878152561,0.997319316611813,1.05050804431090,1.01424299007943,1.03965571749010,1.00313901814010,1.05301959072335,1.01704250703788,1.05809806181616,1.02357814798761,1.31483120932804,1.21111505264918,0.942347056540385,0.925019050140336,0.961782971634532,0.936607318620047,0.960033586814537,0.934807022529764,0.978733617472501,0.955852102428261,0.955911613784748,0.931118738562057,0.970787405473375,0.953775060432683,0.964705605637203,0.945202246004638,0.983836579961716,0.884296826991382,0.955166235803206,0.930175146511089,0.969332029953984,0.950727782794455,0.962403694081432,0.944746749603768,0.981491927667051,0.885571625124068,0.962295388844616,0.941208364100859,0.977588116331144,0.834170222295060,1.00049399187436,0.976793770877953,1.08893287175070,0.141227280443693];

% s = 0.5 * ones(1, length(m));
% start = zeros(10, length(m));
% j = 1;
% for i = 0.5:0.05:1.95
%     start(j,:) = rand(1,length(m));%i*m;
%     j = j+1;
% end

% evaluate the objective fcn at each of these m points
% f=zeros(10,1);
% 
% for jj=1:10
%     xjj_size = size(start(jj,:));
%     input = start(jj,:);
%     f(jj)= utilitycalc(input);
% end
% 
% ind=find(f==min(f));

% %
% gamma=1; % a parameter in RBF (phi function)
% deg=-1; % in the RBF model, p(x) is set to zero.
% % set the mono decreasing lambda sequence
% Ls=[10*0.3.^(0:2),0];
% %
% % set the global minimizer guesstimate
% xstar=2*(2*rand(dim,1)-1);
% %n
% % gradient used?
% useg=1; % gradient is not used
% %
% % set the function used in the RBF approximation
% method='cubic';
% 
% [fmin1,xmin1,iter1,xbar_min1,fbar_min1,fcount]=global_RBF_TRM(myfun,x,f,ind,xstar,Ls,method,deg,gamma,useg);

[fmin,xmin,fcount] = quasi_newton(myfun,m); % run Quasi_Newton local optimizer
