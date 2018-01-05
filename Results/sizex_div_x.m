diffx = zeros(1,186);
for i = 2:187
    diffx(i-1) = norm(xs(i,:) - xs(i-1,:)) / norm(xs(i-1,:));
end