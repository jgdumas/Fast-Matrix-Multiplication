function [e_s,e_1,e_w,e_ks,e_v4,e_vk,e_8,t_s,t_1,t_w,t_ks,t_v4,t_vk,t_8] = error_alternate_real(A,B,n_0)
    scale = comNorm_real(A)*comNorm_real(B);
    E = ec(A,B);
    %E = A*B;
    tic
    C_s = strassen(A,B,n_0);
    t_s = toc;
    e_s = diffe_real(C_s,E)/scale;
    tic
    C_1 = Strassen_alternate(A,B,n_0);
    t_1 = toc;
    e_1 = diffe_real(C_1,E)/scale;
    tic
    C_w = strassenw(A,B,n_0);
    t_w = toc;
    e_w = diffe_real(C_w,E)/scale;
    tic
    C_ks = Winograd_alternate(A,B,n_0);
    t_ks = toc;
    e_ks = diffe_real(C_ks,E)/scale;
    tic
    C_v4 = DPS(A,B,n_0);
    t_v4 = toc;
    e_v4 = diffe_real(C_v4,E)/scale;
    tic
    C_vk = DPS_alternate(A,B,n_0);
    t_vk = toc;
    e_vk = diffe_real(C_vk,E)/scale;
    tic
    C_8 = conventional(A,B,n_0);
    t_8 = toc;
    e_8 = diffe_real(C_8,E)/scale;
end
