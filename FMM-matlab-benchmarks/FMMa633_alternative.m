function C = FMMa633_alternative(A, B, nmin, peeling, level)
%          Computes the product C = A*B, sparsified via alternative basis.
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level

  U = FMMa633_CoBL(A, nmin, peeling, level);   % Left change of basis
  V = FMMa633_CoBR(B, nmin, peeling, level);   % Right change of basis
  W = FMMa633_mul(U, V, nmin, peeling, level); % Sparse multiplication
  C = FMMa633_ICoB(W, nmin, peeling, level);   % Inverse change of basis
end
