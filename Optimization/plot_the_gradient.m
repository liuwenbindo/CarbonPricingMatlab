figure
X = 1:length(utils);
plot(X,utils(X))
xlabel('Times of Iteration')
ylabel('Utility Value')
%legend('AD gradient','FD gradent')
title('Utility Change in Quasi-Newton Optimization')