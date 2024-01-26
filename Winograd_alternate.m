function C = Winograd_alternate(A, B, nmin)
   U = Winograd_CoB(A, nmin);
   V = Winograd_CoB(B, nmin);
   W = Winograd_mul_alternate(U, V, nmin);
   C = Winograd_ICoB(W, nmin);
end
