function [M_1,M_2] = f(A,B,C,D)
    P_1 = A*C;
    P_2 = B*D;
    P_3 = A*D;
    P_4 = B*C;
    M_1 = P_1 - P_2;
    M_2 = P_3 + P_4;
end