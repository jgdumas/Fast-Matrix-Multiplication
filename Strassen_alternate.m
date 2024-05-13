function C = Strassen_alternate(A, B, nmin)
%          Sparsification of Strassen's algorithm
%          via alternative basis.

%          Reference:
%          G. Beniamini, N. Cheng, O. Holtz, E. Karstadt, O. Schwartz;
%          Sparsifying the Operators of Fast Matrix Multiplication Algorithms
%          https://arxiv.org/abs/2008.03759

   U = Strassen_CoB(A, nmin);              % Left change of basis
   V = Strassen_CoB(B, nmin);              % Right change of basis
   W = Strassen_mul_alternate(U, V, nmin); % Sparse multiplication
   C = Strassen_ICoB(W, nmin);             % Inverse change of basis
end
