%   Comparing accuracy of fast 2x2-recursive matrix multiplications
%   Reference:
%   J-G. Dumas, C. Pernet, A. Sedoglavic
%   Strassen's algorithm is not optimally accurate, Feb. 2024
%   ISSAC 2024, Raleigh, NC USA, pp. 254-263.
%   https://hal.science/hal-04441653

rng(1);
t = 49;   % Number of points
num = 9;  % Number of runs for each point
n_0 = 1;  % Basecase of the recursion
n_b = 20; % Smallest matrix dimension
inc = 10; % Matrix dimension increment

E_W = zeros(t,1);E_S = zeros(t,1);E_4 = zeros(t,1);E_J = zeros(t,1);
E_8 = zeros(t,1);E_5 = zeros(t,1);E_3 = zeros(t,1);E_R = zeros(t,1);
T_W = zeros(t,1);T_S = zeros(t,1);T_4 = zeros(t,1);T_J = zeros(t,1);
T_8 = zeros(t,1);T_5 = zeros(t,1);T_3 = zeros(t,1);T_R = zeros(t,1);
N = zeros(t,1);m = zeros(t,1);

parfor l = 1:t
    n = n_b + inc*l;
    m(l) = n;    disp(l);
    for k = 1:num
        A = gen_mat_SW(n);        N(l) = size(A,2);
        B = gen_mat_SW(n);
        [e_w,e_s,e_v4,e_jg,e_ra,e_j5,e_a3,e_8,t_w,t_s,t_v4,t_jg,t_ra,t_j5,t_a3,t_8] = error_2x2_real(A,B,n_0);
        E_W(l) = E_W(l) + e_w/num;        E_4(l) = E_4(l) + e_v4/num;       E_J(l) = E_J(l) + e_jg/num;
        E_S(l) = E_S(l) + e_s/num;        E_8(l) = E_8(l) + e_8/num;        E_5(l) = E_5(l) + e_j5/num;
        E_3(l) = E_3(l) + e_a3/num;       E_R(l) = E_R(l) + e_ra/num;
        T_W(l) = T_W(l) + t_w/num;        T_4(l) = T_4(l) + t_v4/num;       T_J(l) = T_J(l) + t_jg/num;
        T_S(l) = T_S(l) + t_s/num;        T_8(l) = T_8(l) + t_8/num;        T_5(l) = T_5(l) + t_j5/num;
        T_3(l) = T_3(l) + t_a3/num;       T_R(l) = T_R(l) + t_ra/num;
    end
end

save('accuracy_2x2_real','N','m','E_W','E_S','E_J','E_5','E_3','E_4','E_8','E_R','T_W','T_S','T_J','T_5','T_3','T_4','T_8','T_R');
save('accuracy_2x2_real.txt','N','m','E_W','E_S','E_J','E_5','E_3','E_4','E_8','E_R','T_W','T_S','T_J','T_5','T_3','T_4','T_8','T_R',"-ascii","-double");
