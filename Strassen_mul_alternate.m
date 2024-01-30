function C = Strassen_mul_alternate(A, B, nmin)
%          Sparse multiplication of the sparsification
%          of Strassen's algorithm via alternate basis.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
   m = n/2; i = 1:m; j = m+1:n;

   S1 = A(i,i) - A(j,j);
   S2 = A(i,i) + A(j,i);
   S3 = A(i,j) - A(i,i);

   T1 = B(i,i) - B(j,j);
   T2 = B(i,i) + B(j,i);
   T3 = B(i,j) - B(i,i);

   M1 = Strassen_mul_alternate( A(i,i), B(i,i), nmin);
   M2 = Strassen_mul_alternate( A(i,j), T1, nmin);
   M3 = Strassen_mul_alternate( S1, B(j,i), nmin);
   M4 = Strassen_mul_alternate( A(j,j), T3, nmin);
   M5 = Strassen_mul_alternate( S2, B(j,j), nmin);
   M6 = Strassen_mul_alternate( S3, T2, nmin);
   M7 = Strassen_mul_alternate( A(j,i), B(i,j), nmin);

   C11 = M1+M4-M5+M7;
   C12 = M3+M5;
   C21 = M2+M4;
   C22 = M6-M7;

   C = [ C11 C12; C21 C22 ];

end
