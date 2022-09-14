rng(1);
n = 64;
m = 25;
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
        W_1 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        W_2 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        W_3 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        W_4 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        W_5 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        W_6 = gen_mat(n,cond_num)+1j*gen_mat(n,cond_num);
        X = (rand(n,m)-0.5*ones(n,m)) + 1j*(rand(n,m)-0.5*ones(n,m)); %m is the number of test data
        [e_R,e_G,e_N,t_R,t_G,t_N] = error_nn_6(W_1,W_2,W_3,W_4,W_5,W_6,X);
        T_R(l) = T_R(l) + t_R/10;
        E_R(l) = E_R(l) + e_R/10;
        T_G(l) = T_G(l) + t_G/10;
        E_G(l) = E_G(l) + e_G/10;
        T_N(l) = T_N(l) + t_N/10;
        E_N(l) = E_N(l) + e_N/10;
    end
end
save('accuracy_relu_neural_64_depth6','T_R','T_G','T_N','E_R','E_G','E_N');
