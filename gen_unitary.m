function U = gen_unitary(n)
    A = rand(n) + 1j*rand(n);
    [U,R] = qr(A);
end