function [e_f,e_a,e_c,e_d,e_b,e_e,t_f,t_a,t_c,t_d,t_b,t_e] = error_3x6x3_real(A,B,n_0)
	% f fast classical 3x6x3 taken from fmm.univ-lille.fr
	% a fast and accurate
	% c conventional
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_f = FMM_3_6_3(A,B,n_0,1,3);
    t_f = toc;
    e_f = diffe_real(C_f,E)/scale;
    tic
    C_a = FMMa_3_6_3(A,B,n_0,1,3);
    t_a = toc;
    e_a = diffe_real(C_a,E)/scale;
    tic
    C_b = FMM363_alternative(A,B,n_0,1,3);
    t_b = toc;
    e_b = diffe_real(C_b,E)/scale;
    tic
    C_e = FMMa363_alternative(A,B,n_0,1,3);
    t_e = toc;
    e_e = diffe_real(C_e,E)/scale;
    tic
    C_c = recursive(A,B,3,6,3);
    t_c = toc;
    e_c = diffe_real(C_c,E)/scale;
    tic
    C_d = A*B;
    t_d = toc;
    e_d = diffe_real(C_d,E)/scale;
end
