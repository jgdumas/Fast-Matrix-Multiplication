function C = FMMa(A, B, nmin, peeling, level)
%          Computes C = A*B, with fast accurate cores: <6;3;3>,<3;3;6>,<3;6;3>.
%          nmin   : threshold switch to conventional.
%          peeling: static peeling (1, by default) or dynamic (2).
%          level  : logging level.
  if nargin < 3, nmin = 3; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 3; end   % Verbose level
  [m,k] = size(A); [k2,n] = size(B);
  if (k2 ~= k), error('Incompatible matrix dimensions.'); end
  if (m >= max(k,n))
    C=FMMa_6_3_3(A,B,nmin, peeling, level);
  elseif (n >= max(m,k))
    C=FMMa_3_3_6(A,B,nmin, peeling, level);
  else
    C=FMMa_3_6_3(A,B,nmin, peeling, level);
  end
end
