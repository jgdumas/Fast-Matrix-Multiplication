function C = DPS_alternate(A, B, nmin)
   U = DPS_CoBL(A, nmin);
   V = DPS_CoBR(B, nmin);
   W = DPS_mul(U, V, nmin);
   C = DPS_ICoB(W, nmin);
end
