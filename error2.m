function [e_R,e_G,e_N,t_R,t_G,t_N] = error2(A,B,C,D)
    scale = comNorm(A,B)*comNorm(C,D);
    [E_1,E_2] = exact(A,B,C,D);
    tic
    [MR_1,MR_2] = f(A,B,C,D);
    t_R = toc;
    tic
    [MG_1,MG_2] = g(A,B,C,D);
    t_G = toc;
    tic
    [MN_1,MN_2] = h(A,B,C,D);
    t_N = toc;
    e_R = diffe(MR_1,MR_2,E_1,E_2)/scale;
    e_G = diffe(MG_1,MG_2,E_1,E_2)/scale;
    e_N = diffe(MN_1,MN_2,E_1,E_2)/scale;
end
