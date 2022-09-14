rng(1);
t = 50;
T_S = zeros(t,1);
T_W = zeros(t,1);
T_R = zeros(t,1);
E_S = zeros(t,1);
E_W = zeros(t,1);
E_R = zeros(t,1);
n_0 = 1;

for l = 1:t
    disp(l);
    n = 20 + 2*l;
    for k = 1:10
        A = gen_mat_SW(n);
        B = gen_mat_SW(n);
        A_2 = gen_mat_SW(n);
        B_2 = gen_mat_SW(n);
        A = A + 1j*A_2;
        B = B + 1j*B_2;
        [e_S,e_W,e_R,t_S,t_W,t_R] = error_complex_2(A,B,n_0);
        T_S(l) = T_S(l) + t_S/10;
        E_S(l) = E_S(l) + e_S/10;
        T_W(l) = T_W(l) + t_W/10;
        E_W(l) = E_W(l) + e_W/10;
        T_R(l) = T_R(l) + t_R/10;
        E_R(l) = E_R(l) + e_R/10;
    end
end

save('accuracy_SWC_complex','T_S','T_W','T_R','E_S','E_W','E_R');