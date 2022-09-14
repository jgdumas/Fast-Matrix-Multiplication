rng(1);
n = 256;
T_R = zeros(20,1);
T_G = zeros(20,1);
T_N = zeros(20,1);
E_R = zeros(20,1);
E_G = zeros(20,1);
E_N = zeros(20,1);

for l = 1:20
    disp(l);
    cond_num = floor(2^(l+33));
    for k = 1:10
        U = gen_unitary(n);
        A = real(U);
        B = imag(U);
        C = gen_mat(n,cond_num);
        D = gen_mat(n,cond_num);
        [e_R,e_G,e_N,t_R,t_G,t_N] = error2(A,B,C,D);
        T_R(l) = T_R(l) + t_R/10;
        E_R(l) = E_R(l) + e_R/10;
        T_G(l) = T_G(l) + t_G/10;
        E_G(l) = E_G(l) + e_G/10;
        T_N(l) = T_N(l) + t_N/10;
        E_N(l) = E_N(l) + e_N/10;
    end
end
save('accuracy_unitary_256','T_R','T_G','T_N','E_R','E_G','E_N');
