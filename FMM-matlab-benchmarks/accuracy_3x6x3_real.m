%	Comparing accuracy of fast 3x6x3-recursive matrix multiplications
%	Authors: J-G. Dumas, C. Pernet, A. Sedoglavic

rng(1);
t = 4;    % Number of points
num = 9; % Number of runs for each point
n_0 = 1;  % Basecase of the recursion

E_f = zeros(t,1);E_a = zeros(t,1);E_c = zeros(t,1);E_d = zeros(t,1);
T_f = zeros(t,1);T_a = zeros(t,1);T_c = zeros(t,1);T_d = zeros(t,1);
M = zeros(t,1);
N = zeros(t,1);

% A starts with 3x6 and increase (3^l)x(6^l)
% B starts with 6x3 and increase (6^l)x(3^l)

%parfor l = 1:t
for l = 1:t
    disp(l);
    for k = 1:num
	A = gen_mat_rect(3^l,6^l); M(l) = size(A,1);
	B = gen_mat_rect(6^l,3^l);   N(l) = size(B,2);
    disp(M(l));
    disp(N(l));
	[e_f,e_a,e_c,e_d,t_f,t_a,t_c,t_d] = error_3x6x3_real(A,B,n_0);
	E_f(l) = E_f(l) + e_f/num; T_f(l) = T_f(l) + t_f/num;
	E_a(l) = E_a(l) + e_a/num; T_a(l) = T_f(l) + t_f/num;
	E_c(l) = E_c(l) + e_c/num; T_c(l) = T_c(l) + t_c/num;
	E_d(l) = E_d(l) + e_d/num; T_d(l) = T_d(l) + t_d/num;
    end
end

save('accuracy_3x6x3_real','M','N','E_f','E_a','E_c','E_d');
save('accuracy_3x6x3_real.txt','M','N','E_f','E_a','E_c','E_d','T_f','T_a','T_c','T_d',"-ascii","-double");
