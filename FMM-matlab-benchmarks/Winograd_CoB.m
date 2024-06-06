function M = Winograd_CoB(A, nmin)
%          Left & Right change of basis of the sparsification
%          of Winograd's algorithm via alternative basis.

if nargin < 2, nmin = 8; end
n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end
if n <= nmin
   M=A;
else
   m = n/2; i = 1:m; j = m+1:n;

   M1 = Winograd_CoB(A(i,i), nmin);
   M2 = Winograd_CoB(A(i,j), nmin);
   M3 = Winograd_CoB(A(j,i), nmin);
   M4 = Winograd_CoB(A(j,j), nmin);

   S1 = M1;
   S3 = M4-M3;
   S2 = M2+S3;
   S4 = M2+M4;

   M = [ S1 S2; S3 S4];
end
