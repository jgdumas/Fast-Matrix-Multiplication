function [e_R,e_G,e_N,t_R,t_G,t_N] = error_nn_6(W_1,W_2,W_3,W_4,W_5,W_6,X)
    [E_1,E_2] = neural_network_6(W_1,W_2,W_3,W_4,W_5,W_6,X,@ec);
    scale = comNorm(double(E_1),double(E_2));
    tic
    [MR_1,MR_2] = neural_network_6(W_1,W_2,W_3,W_4,W_5,W_6,X,@fc);
    t_R = toc;
    tic
    [MG_1,MG_2] = neural_network_6(W_1,W_2,W_3,W_4,W_5,W_6,X,@gc);
    t_G = toc;
    tic
    [MN_1,MN_2] = neural_network_6(W_1,W_2,W_3,W_4,W_5,W_6,X,@hc);
    t_N = toc;
    e_R = diffe(MR_1,MR_2,E_1,E_2)/scale;
    e_G = diffe(MG_1,MG_2,E_1,E_2)/scale;
    e_N = diffe(MN_1,MN_2,E_1,E_2)/scale;
end