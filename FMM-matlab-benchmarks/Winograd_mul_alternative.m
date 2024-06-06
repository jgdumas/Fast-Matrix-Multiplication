function C = Winograd_mul_alternative(A, B, nmin)
%          Sparse multiplication of the sparsification
%          of Winograd's algorithm via alternative basis.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
   m = n/2; i = 1:m; j = m+1:n;

   S1 = A(i,j) - A(j,i);
   S2 = A(i,j) - A(i,i);
   S3 = A(j,j) - A(i,j);

   T1 = B(j,j) - B(i,j);
   T2 = B(i,j) - B(j,i);
   T3 = B(i,j) - B(i,i);

   M1 = Winograd_mul_alternative( A(j,j), B(j,j), nmin);
   M2 = Winograd_mul_alternative( A(j,i), B(j,i), nmin);
   M3 = Winograd_mul_alternative( A(i,j), B(i,j), nmin);
   M4 = Winograd_mul_alternative( A(i,i), B(i,i), nmin);
   M5 = Winograd_mul_alternative( S1, T1, nmin);
   M6 = Winograd_mul_alternative( S2, T2, nmin);
   M7 = Winograd_mul_alternative( S3, T3, nmin);

   C11 = M4+M5;
   C12 = M3+M5-M6+M7;
   C21 = M2+M7;
   C22 = M1-M6;

   C = [ C11 C12; C21 C22 ];

end
