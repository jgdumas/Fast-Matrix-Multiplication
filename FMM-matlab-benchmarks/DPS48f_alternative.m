function C = DPS48f_alternative(A, B, nmin, peeling, level)
%          Computes the product C = A*B, sparsified via alternative basis.
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level

  U = DPS48f_CoBL(A, nmin, peeling, level);   % Left change of basis
  V = DPS48f_CoBR(B, nmin, peeling, level);   % Right change of basis
  W = DPS48f_mul(U, V, nmin, peeling, level); % Sparse multiplication
  C = DPS48f_ICoB(W, nmin, peeling, level);   % Inverse change of basis
end
