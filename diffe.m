function d = diffe(M_1,M_2,E_1,E_2)
    d = comNorm(M_1 - double(E_1), M_2 - double(E_2));
end