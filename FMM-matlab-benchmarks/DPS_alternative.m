function C = DPS_alternative(A, B, nmin)
%          Sparsification of DPS's algorithm via alternative basis.

   U = DPS_CoBL(A, nmin);	% Left change of basis
   V = DPS_CoBR(B, nmin);	% Right change of basis
   W = DPS_mul(U, V, nmin); % Sparse multiplication
   C = DPS_ICoB(W, nmin);   % Inverse change of basis
end
