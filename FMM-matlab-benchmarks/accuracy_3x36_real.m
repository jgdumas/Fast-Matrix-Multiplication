%   Comparing accuracy of fast 2x2-recursive matrix multiplications
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

E_f = zeros(t,1);E_fa = zeros(t,1);E_c = zeros(t,1);
T_f = zeros(t,1);T_fa = zeros(t,1);T_c = zeros(t,1);
N = zeros(t,1);m = zeros(t,1);

parfor l = 1:t
    n = n_b + inc*l;
    m(l) = n;    disp(l);
    for k = 1:num
        A = gen_mat_SW(n);        N(l) = size(A,2);  % should this be modified ? 
        B = gen_mat_SW(n);
        [e_f,e_fa,e_c,t_f,t_fa,t_c] = error_3x3x6_real(A,B,n_0);
        %E_f(l) = E_f(l) + e_f/num;      T_f(l) = T_f(l) + t_f/num;
        %E_fa(l) = E_fa(l) + e_fa/num; T_fa(l) = T_fa(l) + t_fa/num;      
        E_c(l) = E_c(l) + e_c/num; T_c(l) = T_c(l) + t_c/num;  
    end
end

save('accuracy_2x2_real','N','m','E_f','E_fa','E_c');
save('accuracy_2x2_real.txt','N','m','E_f','E_fa','E_c''T_f','T_fa','T_c',"-ascii","-double");
