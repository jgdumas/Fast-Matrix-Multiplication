function C = conventional_3x3x6(A, B, nmin)
%Conventional	Conventional matrix multiplication algorithm.
%          	C = conventional_3x3x6(A, B, NMIN), where A and B are matrices of dimension
%          	a power of ?, computes the product C = A*B.
%          	conventional algorithm is used recursively until dimension <= NMIN
%          	is reached, at which point standard multiplication is used.
%          	The default is NMIN = 6 (which minimizes the total number of
%          	operations).

if nargin < 3, nmin = 6; end

ca = size(A);cb = size(B);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if ca(2) <= nmin
   C = A*B;
else
   i = 1:(ca/3); j = (ca/3)+1:2*(ca/3); k:=2*(ca/3)+1:ca;
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
