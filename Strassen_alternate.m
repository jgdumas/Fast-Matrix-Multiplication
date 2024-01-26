function C = Strassen_alternate(A, B, nmin)
   U = Strassen_CoB(A, nmin);
   V = Strassen_CoB(B, nmin);
   W = Strassen_mul_alternate(U, V, nmin);
   C = Strassen_ICoB(W, nmin);
end
