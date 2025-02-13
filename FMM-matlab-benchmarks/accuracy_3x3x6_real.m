%	Comparing accuracy of fast 3x3x6-recursive matrix multiplications
%	Reference:
%	J-G. Dumas, C. Pernet, A. Sedoglavic
%	Next level

rng(1);
t = 3;   % Number of points
num = 5;  % Number of runs for each point
n_0 = 1;  % Basecase of the recursion

E_f = zeros(t,1);E_a = zeros(t,1);E_c = zeros(t,1);
T_f = zeros(t,1);T_a = zeros(t,1);T_c = zeros(t,1);
M = zeros(t,1);
N = zeros(t,1);

% A starts with 3x3 and increase (3^l)x(3^l)
% B starts with 3x6 and increase (3^l)x(6^l)

%parfor l = 1:t
for l = 1:t
    disp(l);
    for k = 1:num
	A = gen_mat_rect(3^l,3^l); M(l) = size(A,1);
	B = gen_mat_rect(3^l,6^l);   N(l) = size(B,2);
    disp(M(l));
    disp(N(l));
	[e_f,e_a,e_c,t_f,t_a,t_c] = error_3x3x6_real(A,B,n_0);
	E_f(l) = E_f(l) + e_f/num; T_f(l) = T_f(l) + t_f/num;
	E_a(l) = E_a(l) + e_a/num; T_a(l) = T_f(l) + t_f/num;
	E_c(l) = E_c(l) + e_c/num; T_c(l) = T_c(l) + t_c/num;
    end
end

save('accuracy_3x3x6_real','M','N','E_f','E_a','E_c');
save('accuracy_3x3x6_real.txt','M','N','E_f','E_a','E_c','T_f','T_a','T_c',"-ascii","-double");
