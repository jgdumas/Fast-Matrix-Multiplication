function [e_f,e_a,e_c,t_f,t_a,t_c] = error_3x3x6_real(A,B,n_0)
	% f fast classical 3x3x6 taken from fmm.univ-lille.fr
	% a fast and accurate
	% c conventional
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_f = FMM_3_3_6(A,B,n_0,1);
    t_f = toc;
    e_f = diffe_real(C_f,E)/scale;
    tic
    C_a = FMMa_3_3_6(A,B,n_0,1);
    t_a = toc;
    e_a = diffe_real(C_a,E)/scale;
    tic
    C_c = conventional_3x3x6(A,B,n_0);
    t_c = toc;
    e_c = diffe_real(C_c,E)/scale;
end
