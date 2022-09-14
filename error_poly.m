function [e_R,e_G,e_N,t_R,t_G,t_N] = error_poly(X,b)
    [E_1,E_2] = poly_eval(sym(X),sym(b),@ec);
    scale = comNorm(double(E_1),double(E_2));
    tic
    [MR_1,MR_2] = poly_eval(X,b,@fc);
    t_R = toc;
    tic
    [MG_1,MG_2] = poly_eval(X,b,@gc);
    t_G = toc;
    tic
    [MN_1,MN_2] = poly_eval(X,b,@hc);
    t_N = toc;
    e_R = diffe(MR_1,MR_2,E_1,E_2)/scale;
    e_G = diffe(MG_1,MG_2,E_1,E_2)/scale;
    e_N = diffe(MN_1,MN_2,E_1,E_2)/scale;
end