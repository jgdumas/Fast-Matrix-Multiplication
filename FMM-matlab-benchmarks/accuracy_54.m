%	Comparing accuracy of fast 54-recursive matrix multiplications
%	Authors: J-G. Dumas, C. Pernet, A. Sedoglavic

rng(1);
t = 5;   % Number of points
num = 9; % Number of runs for each point
n_0 = 1; % Basecase of the recursion

E_f = zeros(t,1);E_a = zeros(t,1);E_c = zeros(t,1);E_d = zeros(t,1);
T_f = zeros(t,1);T_a = zeros(t,1);T_c = zeros(t,1);T_d = zeros(t,1);
M = zeros(t,1); K = zeros(t,1); N = zeros(t,1);

% A starts with 6x3 and B starts with 3x3


%parfor l = 1:t
for l = 1:t
  disp(l);
  for k = 1:num
    q=floor(l/3);
    mu=2^q*3^l; ku=2^q*3^l; nu=2^q*3^l;
    if (rem(l,3)>=1), mu=mu*2; end
    if (rem(l,3)>=2), nu=nu*2; end
	A = gen_mat_rect(mu,ku); M(l) = size(A,1); K(l) = size(A,2);
	B = gen_mat_rect(ku,nu); N(l) = size(B,2);
    disp([M(l),K(l),N(l)]);
	[e_f,e_a,e_c,e_d,t_f,t_a,t_c,t_d] = error_FMM(A,B,n_0,2);
	E_f(l) = E_f(l) + e_f/num; T_f(l) = T_f(l) + t_f/num;
	E_a(l) = E_a(l) + e_a/num; T_a(l) = T_f(l) + t_f/num;
	E_c(l) = E_c(l) + e_c/num; T_c(l) = T_c(l) + t_c/num;
	E_d(l) = E_d(l) + e_d/num; T_d(l) = T_d(l) + t_d/num;
  end
end

save('accuracy_54','M','N','E_f','E_a','E_c','E_d');
save('accuracy_54.txt','M','N','E_f','E_a','E_c','E_d','T_f','T_a','T_c','T_d',"-ascii","-double");
