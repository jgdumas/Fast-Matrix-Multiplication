function C = Winograd_alternate(A, B, nmin)
%          Sparsification of Winograd's algorithm
%          via alternate basis.

%          Reference:
%          E. Karstadt, O. Schwartz; Matrix multiplication, a little faster.
%          SPAA 2017, pages 101--110, https://doi.org/10.1145/3087556.3087579

   U = Winograd_CoB(A, nmin);              % Left change of basis
   V = Winograd_CoB(B, nmin);              % Right change of basis
   W = Winograd_mul_alternate(U, V, nmin); % Sparse multiplication
   C = Winograd_ICoB(W, nmin);             % Inverse change of basis
end
