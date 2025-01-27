%   	Comparing accuracy of fast 3x3x6-recursive matrix multiplications
%   	Reference:
%   	J-G. Dumas, C. Pernet, A. Sedoglavic
%	Next level

rng(1);
t = 4;   % Number of points
num = 9;  % Number of runs for each point
n_0 = 3;  % Basecase of the recursion

E_f = zeros(t,1);E_fa = zeros(t,1);E_c = zeros(t,1);
T_f = zeros(t,1);T_fa = zeros(t,1);T_c = zeros(t,1);
N = zeros(t,1);m = zeros(t,1);

% A starts with 3x3 and increase (3+l)x(3+l)
% B starts with 3x6 and increase (3+l)x(3+2*l)
n = 6; % Smallest first matrix dimension
inc = 6; % Matrix dimension increment

%parfor l = 1:t
for l = 1:t
    n = n + inc*l;
    m(l) = n_1;    disp(l);
    for k = 1:num
        A = gen_mat_SW(n);        N(l) = size(A,2);  
        B = gen_mat_SW(n);
        [e_c,t_c] = error_3x3x6_real(A,B,n_0);
        [e_f,e_c,t_f,t_c] = error_3x3x6_real(A,B,n_0);
        %[e_f,e_fa,e_c,t_f,t_fa,t_c] = error_3x3x6_real(A,B,n_0);
        E_f(l) = E_f(l) + e_f/num;      T_f(l) = T_f(l) + t_f/num;
        %E_fa(l) = E_fa(l) + e_fa/num; T_fa(l) = T_fa(l) + t_fa/num;      
        E_c(l) = E_c(l) + e_c/num; T_c(l) = T_c(l) + t_c/num;  
    end
end

save('accuracy_3x3x6_real','N','m','E_f','E_fa','E_c');
save('accuracy_3x3x6_real.txt','N','m','E_f','E_fa','E_c','T_f','T_fa','T_c',"-ascii","-double");
