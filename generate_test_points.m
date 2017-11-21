x = 5 * ones(1,63);
% for i = 1:63
% if mod(i,2) == 0
% x(i) = x(i) +0.1;
% end
% end
for i = 2:63

x(i) = x(i-1) -0.05;

end
utilitycalc(x)