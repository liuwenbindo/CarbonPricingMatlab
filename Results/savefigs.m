for k = 1:10
    figure
    fig = bar(eig(H_random(:,:,k)));
    name = strcat(num2str(k),'.jpg');
    saveas(fig,name);
end