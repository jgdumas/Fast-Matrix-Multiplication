function M = Strassen_ICoB(A, nmin)
%          Inverse change of basis of the sparsification
%          of Strassen's algorithm via alternative basis.

if nargin < 2, nmin = 8; end
n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end
if n <= nmin
   M=A;
else
   m = n/2; i = 1:m; j = m+1:n;

   M1 = Strassen_ICoB(A(i,i), nmin);
   M2 = Strassen_ICoB(A(i,j), nmin);
   M3 = Strassen_ICoB(A(j,i), nmin);
   M4 = Strassen_ICoB(A(j,j), nmin);

   S1 = M1;
   S2 = M2;
   S3 = M3;
   S4 = M1+M2-M3+M4;

   M = [ S1 S2; S3 S4];
end
