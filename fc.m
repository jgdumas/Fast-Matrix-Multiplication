function Z = fc(X,Y)
    [Z_r,Z_c] = f(real(X),imag(X),real(Y),imag(Y));
    Z = Z_r + 1j*Z_c;
end