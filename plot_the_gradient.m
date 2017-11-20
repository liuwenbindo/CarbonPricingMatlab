figure
X = 1:63;
plot(X,g1,X,g2)
xlabel('Dimension')
ylabel('Gradient Value')
legend('AD gradient','FD gradent')
title('Comparison of AD/FD Gradient at the New Optimal Point m*')