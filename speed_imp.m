rng(1);
T_R = zeros(50,1);
T_G = zeros(50,1);
T_N = zeros(50,1);
for l = 1:50
    disp(l);
    n = 2000 + 100*l;
    [t_R,t_G,t_N] = speed(n);
    T_R(l) = t_R;
    T_G(l) = t_G;
    T_N(l) = t_N;
end
save('speed_general_2100_7000','T_R','T_G','T_N')