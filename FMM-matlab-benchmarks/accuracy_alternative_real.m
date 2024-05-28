%   Comparing accuracy of sparsified 2x2-recursive matrix multiplications
%   Reference:
%   J-G. Dumas, C. Pernet, A. Sedoglavic
%   Strassen's algorithm is not optimally accurate, Feb. 2024
%   https://hal.science/hal-04441653

rng(1);
t = 49;   % Number of points
num = 9;  % Number of runs for each point
n_0 = 1;  % Basecase of the recursion
n_b = 20; % Smallest matrix dimension
inc = 10; % Matrix dimension increment

E_S = zeros(t,1);E_1 = zeros(t,1);E_W = zeros(t,1);E_K = zeros(t,1);
E_4 = zeros(t,1);E_O = zeros(t,1);E_8 = zeros(t,1);
T_S = zeros(t,1);T_1 = zeros(t,1);T_W = zeros(t,1);T_K = zeros(t,1);
T_4 = zeros(t,1);T_O = zeros(t,1);T_8 = zeros(t,1);
N = zeros(t,1);m = zeros(t,1);

parfor l = 1:t
    n = n_b + inc*l;
    m(l) = n; disp(l);
    for k = 1:num
        A = gen_mat_SW(n);        N(l) = size(A,2);
        B = gen_mat_SW(n);
        [e_s,e_1,e_w,e_ks,e_v4,e_vk,e_8,t_s,t_1,t_w,t_ks,t_v4,t_vk,t_8] = error_alternative_real(A,B,n_0);
        E_S(l) = E_S(l) + e_s/num;        E_1(l) = E_1(l) + e_1/num;        E_W(l) = E_W(l) + e_w/num;
        E_4(l) = E_4(l) + e_v4/num;       E_O(l) = E_O(l) + e_vk/num;       E_K(l) = E_K(l) + e_ks/num;
        E_8(l) = E_8(l) + e_8/num;
        T_S(l) = T_S(l) + t_s/num;        T_1(l) = T_1(l) + t_1/num;        T_W(l) = T_W(l) + t_w/num;
        T_4(l) = T_4(l) + t_v4/num;       T_O(l) = T_O(l) + t_vk/num;       T_K(l) = T_K(l) + t_ks/num;
        T_8(l) = T_8(l) + t_8/num;
    end
end

save('accuracy_alternative_real','N','m','E_S','E_1','E_W','E_K','E_4','E_O','E_8','T_S','T_1','T_W','T_K','T_4','T_O','T_8');
save('accuracy_alternative_real.txt','N','m','E_S','E_1','E_W','E_K','E_4','E_O','E_8','T_S','T_1','T_W','T_K','T_4','T_O','T_8',"-ascii","-double");
