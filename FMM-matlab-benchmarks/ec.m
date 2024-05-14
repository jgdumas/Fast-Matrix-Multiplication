function Z = ec(X,Y)
    %[Z_r,Z_c] = exact(real(X),imag(X),real(Y),imag(Y));
    %Z = Z_r + 1j*Z_c;
    Z = sym(X)*sym(Y);
end