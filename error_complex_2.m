function [e_S,e_W,e_R,t_S,t_W,t_R] = error_complex_2(A,B,n_0)
    scale = comNorm(real(A),imag(A))*comNorm(real(B),imag(B));
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
    e_S = diffe(real(C_S),imag(C_S),real(E),imag(E))/scale;
    e_W = diffe(real(C_W),imag(C_W),real(E),imag(E))/scale;
    e_R = diffe(real(C_R),imag(C_R),real(E),imag(E))/scale;
end