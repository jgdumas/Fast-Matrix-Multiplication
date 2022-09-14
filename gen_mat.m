function T = gen_mat(n,cond_num)
    S = randi([1,cond_num - 1],[n-2,1]);
    S = [S;cond_num;1];
    S = diag(S);
    T = hadamard(n)*S*hadamard(n)';
    T = T/n;
end