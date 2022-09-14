function [M_1,M_2] = exact(A,B,C,D)
    [M_1,M_2] = f(sym(A),sym(B),sym(C),sym(D));
end