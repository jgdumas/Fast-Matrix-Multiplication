function [B, R] = gen_mat_dvs(n, S)
    R = S \ randn(n);
    m = 2^(ceil(log2(n)));
    B = [R,zeros(n,m-n);zeros(m-n,m)];
end
