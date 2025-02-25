function C = DPS(A, B, nmin, peeling, level)
%          Computes the product C = A*B, with fast core: <2;2;2>.
%          nmin   : threshold switch to conventional.
%          peeling: static peeling (1, by default) or dynamic (2).
%          level  : logging level.

%          Used recursively until dimension <= NMIN, is reached,
%          at which point standard multiplication is used.
%          Growth-factor: 2\sqrt(2)+16\sqrt(3)\approx 12.06603145.
%          Frobenius norms: \sqrt(10),\sqrt(10),\sqrt(10),
%          7+7+10 additions, 4+4+4 multiplications/divisions.

%          Reference:
%          J-G. Dumas, C. Pernet, A. Sedoglavic
%          Strassen's algorithm is not optimally accurate, Feb. 2024
%          ISSAC 2024, Raleigh, NC USA, pp. 254-263.
%          https://hal.science/hal-04441653

  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  [m,k] = size(A); [k2,n] = size(B);
  if (k2 ~= k), error('Incompatible matrix dimensions.'); end
  if nargin < 5                   % Min level for verbose output
    level = min([floor(log(m/nmin)/log(2)),floor(log(k/nmin)/log(2)),floor(log(n/nmin)/log(2))]);
  end
  C=DPS_2_2_2(A,B,nmin,peeling,level);
end
