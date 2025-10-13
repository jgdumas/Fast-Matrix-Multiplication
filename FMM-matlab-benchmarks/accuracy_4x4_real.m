%   Comparing accuracy of fast 4x4-recursive matrix multiplications
%   Reference:
%   J-G. Dumas, C. Pernet, A. Sedoglavic
%   Strassen's algorithm is not optimally accurate, Feb. 2024
%   ISSAC 2024, Raleigh, NC USA, pp. 254-263.
%   https://hal.science/hal-04441653

rng(1);
t = 60;   % Number of points
num = 7;  % Number of runs for each point
n_0 = 1;  % Basecase of the recursion
n_b = 24; % Smallest matrix dimension
inc = 8; % Matrix dimension increment

E_W = zeros(t,1);E_S = zeros(t,1);E_4 = zeros(t,1);E_B = zeros(t,1);
E_8 = zeros(t,1);E_O = zeros(t,1);E_A = zeros(t,1);E_F = zeros(t,1);
T_W = zeros(t,1);T_S = zeros(t,1);T_4 = zeros(t,1);T_B = zeros(t,1);
T_8 = zeros(t,1);T_O = zeros(t,1);T_A = zeros(t,1);T_F = zeros(t,1);
N = zeros(t,1);m = zeros(t,1);

parfor l = 1:t
    n = n_b + inc*l;
    m(l) = n;    disp(l); disp(n);
    for k = 1:num
        A = gen_mat_SW(n);        N(l) = size(A,2);
        B = gen_mat_SW(n);
        [e_w,e_s,e_v4,e_b,e_o,e_a,e_f,e_8,t_w,t_s,t_v4,t_b,t_o,t_a,t_f,t_8] = error_4x4_real(A,B,n_0);
        E_W(l) = E_W(l) + e_w/num;        E_4(l) = E_4(l) + e_v4/num;
        E_S(l) = E_S(l) + e_s/num;        E_8(l) = E_8(l) + e_8/num; E_B(l) = E_B(l) + e_b/num;
        E_O(l) = E_O(l) + e_o/num; E_A(l) = E_A(l) + e_a/num; E_F(l) = E_F(l) + e_f/num;
        T_W(l) = T_W(l) + t_w/num;        T_4(l) = T_4(l) + t_v4/num;
        T_S(l) = T_S(l) + t_s/num;        T_8(l) = T_8(l) + t_8/num; T_B(l) = T_B(l) + t_b/num;
        T_O(l) = T_O(l) + t_o/num;
        T_A(l) = T_A(l) + t_a/num;
        T_F(l) = T_F(l) + t_f/num;
    end
end
fprintf("Saving ... ");

save('accuracy_4x4_real','N','m','E_W','E_S','E_B','E_O','E_A','E_F','E_4','E_8','T_W','T_S','T_B','T_O','T_A','T_F','T_4','T_8');
save('accuracy_4x4_real.txt','N','m','E_W','E_S','E_B','E_O','E_A','E_F','E_4','E_8','T_W','T_S','T_B','T_O','T_A','T_F','T_4','T_8',"-ascii","-double");
fprintf(" ... done.\n");
