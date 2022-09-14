function [e_S,e_W,e_R,t_S,t_W,t_R] = error_real_2(A,B,n_0)
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    tic
    C_S = strassen(A,B,n_0);
    t_S = toc;
    tic
    C_W = strassenw(A,B,n_0);
    t_W = toc;
    tic
    C_R = conventional(A,B,n_0);
    t_R = toc;
    e_S = diffe_real(C_S,E)/scale;
    e_W = diffe_real(C_W,E)/scale;
    e_R = diffe_real(C_R,E)/scale;
end