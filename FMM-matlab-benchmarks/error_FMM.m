function [e_f,e_a,e_c,e_d,t_f,t_a,t_c,t_d] = error_FMM(A,B,n_0,peeling)
	% f fast from fmm.univ-lille.fr
	% a fast and accurate
	% c conventional
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_f = FMM(A,B,n_0,peeling,3);
    t_f = toc;
    e_f = diffe_real(C_f,E)/scale;
    tic
    C_a = FMMa(A,B,n_0,peeling,3);
    t_a = toc;
    e_a = diffe_real(C_a,E)/scale;
    tic
    C_c = recursive(A,B,128,128,128);
    t_c = toc;
    e_c = diffe_real(C_c,E)/scale;
    tic
    C_d = A*B;
    t_d = toc;
    e_d = diffe_real(C_d,E)/scale;
end
