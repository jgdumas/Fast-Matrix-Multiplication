function C = DPS_integral(A, B, nmin, peeling, level)
%          Computes the product C = A*B, with fast core: <2;2;2>.
%          nmin   : threshold switch to conventional.
%          peeling: static peeling (1, by default) or dynamic (2).
%          level  : logging level.
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  [m,k] = size(A); [k2,n] = size(B);
  if (k2 ~= k), error('Incompatible matrix dimensions.'); end
  if nargin < 5                   % Min level for verbose output
    level = min([floor(log(m/nmin)/log(2)),floor(log(k/nmin)/log(2)),floor(log(n/nmin)/log(2))]);
  end
  C=DPS_integral_2_2_2(A,B,nmin,peeling,level);
end
