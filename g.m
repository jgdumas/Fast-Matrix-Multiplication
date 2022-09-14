function [M_1,M_2] = g(A,B,C,D)
    P_1 = A*C;
    P_2 = B*D;
    S_3 = A + B;
    S_4 = C + D;
    P_3 = S_3*S_4;
    M_1 = P_1 - P_2;
    M_2 = P_3 - P_1 - P_2;
end