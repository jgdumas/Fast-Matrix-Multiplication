function [e_w,e_s,e_v4,e_b,e_o,e_a,e_f,e_8,t_w,t_s,t_v4,t_b,t_o,t_a,t_f,t_8] = error_4x4_real(A,B,n_0)
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    %E=A*B;
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
    C_o = DPS48o(A,B,n_0);
    t_o = toc;
    e_o = diffe_real(C_o,E)/scale;
    tic
    C_b = DPS48(A,B,n_0);
    t_b = toc;
    e_b = diffe_real(C_b,E)/scale;
    tic
    C_a = DPS48a_alternative(A,B,n_0,0);
    t_a = toc;
    e_a = diffe_real(C_a,E)/scale;
    tic
    C_f = DPS48f_alternative(A,B,n_0,0);
    t_f = toc;
    e_f = diffe_real(C_f,E)/scale;
    tic
    C_8 = conventional(A,B,n_0);
    t_8 = toc;
    e_8 = diffe_real(C_8,E)/scale;
end
