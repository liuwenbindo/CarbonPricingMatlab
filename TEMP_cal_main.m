function cal_main()
    cal();
end

function r = calx(x1,x2)
    r = x1 * x2;
end

function returnd = cal()
    a = 10;
    b = 20;
    job = batch(@calx, 1, {10,20}, 'Pool',2);
    wait(job);
    load(job);
    result = fetchOutputs(job); 
    returnd = result{1};
end