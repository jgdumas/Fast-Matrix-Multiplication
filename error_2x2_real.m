function [e_w,e_s,e_v4,e_jg,e_ra,e_j5,e_a3,e_8,t_w,t_s,t_v4,t_jg,t_ra,t_j5,t_a3,t_8] = error_2x2_real(A,B,n_0)
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_w = strassenw(A,B,n_0);
    t_w = toc;
    e_w = diffe_real(C_w,E)/scale;
    tic
    C_s = strassen(A,B,n_0);
    t_s = toc;
    e_s = diffe_real(C_s,E)/scale;
    tic
    C_v4 = DPS(A,B,n_0);
    t_v4 = toc;
    e_v4 = diffe_real(C_v4,E)/scale;
    tic
    C_jg = DPS_evenpow(A,B,n_0);
    t_jg = toc;
    e_jg = diffe_real(C_jg,E)/scale;
    tic
    C_ra = DPS_smallrat(A,B,n_0);
    t_ra = toc;
    e_ra = diffe_real(C_ra,E)/scale;
    tic
    C_j5 = DPS_intermediate(A,B,n_0);
    t_j5 = toc;
    e_j5 = diffe_real(C_j5,E)/scale;
    tic
    C_a3 = DPS_integral(A,B,n_0);
    t_a3 = toc;
    e_a3 = diffe_real(C_a3,E)/scale;
    tic
    C_8 = conventional(A,B,n_0);
    t_8 = toc;
    e_8 = diffe_real(C_8,E)/scale;
end
