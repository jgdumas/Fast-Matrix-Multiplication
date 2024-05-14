function C = conventional(A, B, nmin)
%STRASSEN  Strassen's fast matrix multiplication algorithm.
%          C = STRASSEN(A, B, NMIN), where A and B are matrices of dimension
%          a power of 2, computes the product C = A*B.
%          Strassen's algorithm is used recursively until dimension <= NMIN
%          is reached, at which point standard multiplication is used.
%          The default is NMIN = 8 (which minimizes the total number of
%          operations).

%          Reference:
%          V. Strassen, Gaussian elimination is not optimal,
%          Numer. Math., 13 (1969), pp. 354-356.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
   m = n/2; i = 1:m; j = m+1:n;
   P1 = conventional(A(i,i),B(i,i),nmin);
   P2 = conventional(A(i,j),B(j,i),nmin);
   P3 = conventional(A(i,i),B(i,j),nmin);
   P4 = conventional(A(i,j),B(j,j),nmin);
   P5 = conventional(A(j,i),B(i,i),nmin);
   P6 = conventional(A(j,j),B(j,i),nmin);
   P7 = conventional(A(j,i),B(i,j),nmin);
   P8 = conventional(A(j,j),B(j,j),nmin);
   C = [ P1+P2  P3+P4;  P5+P6  P7+P8 ];
end