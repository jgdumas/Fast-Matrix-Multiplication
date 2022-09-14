function [M_1,M_2] = h(A,B,C,D)
    r = sqrt(3);
    r_1 = 1/r;
    r_2 = r/2;
    B_1 = r_1 * B;
    D_1 = r_1 * D;
    S_1 = A + B_1;
    S_2 = C + D_1;
    S_3 = A - B_1;
    S_4 = C - D_1;
    P_1 = S_1*S_2;
    P_2 = S_3*S_4;
    P_3 = B*D;
    M_1 = 1/2 * (P_1 + P_2 - 8/3 * P_3);
    M_2 = r_2 * (P_1 - P_2);
end