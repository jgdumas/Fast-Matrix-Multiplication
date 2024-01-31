function C = DPS_mul(A, B, nmin)
%          Sparse multiplication of the sparsification
%          of DPS's algorithm via alternate basis.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
   m = n/2; i = 1:m; j = m+1:n;

   S1 = A(j,i) - A(j,j);
   S2 = A(i,i) + A(j,j);
   S3 = A(i,j) + A(j,j);
   T1 = B(j,i) - B(j,j);
   T2 = B(i,i) - B(j,j);
   T3 = B(i,j) + B(j,j);

   M1 = DPS_mul( S1, B(i,i), nmin);
   M2 = DPS_mul( A(j,i), -B(i,j), nmin);
   M3 = DPS_mul( A(i,j), B(j,i), nmin);
   M4 = DPS_mul( -A(i,i), T1, nmin);
   M5 = DPS_mul( A(j,j), B(j,j), nmin);
   M6 = DPS_mul( S2, T2, nmin);
   M7 = DPS_mul( S3, T3, nmin);

   C11 = M6+M7;
   C12 = M3-M1;
   C21 = M2+M4;
   C22 = M1+M4+M5+M6;

   C = [ C11 C12; C21 C22 ];

end
