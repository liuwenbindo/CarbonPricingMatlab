condH= zeros(1,10);
for i=1:10
    condH(i) = cond(H_random(:,:,i));
end