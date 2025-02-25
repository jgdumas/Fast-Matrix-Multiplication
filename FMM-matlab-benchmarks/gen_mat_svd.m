function [A,T] = gen_mat_svd(c,n)
    T = gallery('randsvd', n, c, 3, n, n, 1);
    m = 2^(ceil(log2(n)));
    A = [T,zeros(n,m-n);zeros(m-n,m)];
end
