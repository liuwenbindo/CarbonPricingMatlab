g = zeros(1, 187);
myfun = ADfun('utilitycalc',1);
for i = 1:187
   [~,gg] = feval(myfun, xs(i,:));
   g(i) = norm(gg);
   gg = 0;
end