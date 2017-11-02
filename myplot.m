figure
X = 1:63;
plot(X,G1(5,:),X,G2(5,:))
xlabel('Dimension')
ylabel('Gradient Value')
legend('AD gradient','FD gradent')
title('Comparison of AD/FD Gradient in [0.9, 0.7, 0.9, 0.7...]')