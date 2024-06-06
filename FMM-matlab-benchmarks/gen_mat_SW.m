function T = gen_mat_SW(n)
    T = 2*rand(n)-1;
    %T = randn(n); use this line if normal distribution is used
    m = 2^(ceil(log2(n)));
    T = [T,zeros(n,m-n);zeros(m-n,m)];
end