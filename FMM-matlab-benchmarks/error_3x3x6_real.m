function [e_f,e_c,t_f,t_c] = error_3x3x6_real(A,B,n_0)
%function [e_f,e_fa,e_c,t_f,t_fa,t_c] = error_3x3x6_real(A,B,n_0)
	% f fast classical 3x3x6 taken from fmm.univ-lille.fr
	% fa fast and accurate (ours)
	% c conventional
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_f = 3x3x6_f(A,B,n_0);
    t_f = toc;
    e_f = diffe_real(C_f,E)/scale;
    tic
    %C_fa = 3x3x6_fa(A,B,n_0);
    %t_fa = toc;
    %e_fa = diffe_real(C_fa,E)/scale;
    %tic
    C_c = conventional_3x3x6(A,B,n_0);
    t_c = toc;
    e_c = diffe_real(C_c,E)/scale;
end
