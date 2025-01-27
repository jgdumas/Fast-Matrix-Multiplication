function C = conventional_3x3x6(A, B, nmin)
%Conventional	Conventional matrix multiplication algorithm.
%          	C = conventional_3x3x6(A, B, NMIN), where A and B are matrices of dimension
%          	a power of ?, computes the product C = A*B.
%          	conventional algorithm is used recursively until dimension <= NMIN
%          	is reached, at which point standard multiplication is used.
%          	The default is NMIN = 3 (which minimizes the total number of
%          	operations).

if nargin < 3, nmin = 3; end

ca = size(A,1);cb = size(B,2);
if ((ca <= nmin) | (cb <=nmin))
   C = A*B;
else
if ((rem(ca,3) ~= 0) | (rem(size(A,2),3) ~=0))
   error('The first matrix dimension must be a multiple of 3.')
end
if ((rem(cb,3) ~= 0) | (rem(size(B,2),6) ~=0))
   error('The second matrix dimension must be a multiple of [3,6].')
end

   atom = ca/3;
   i = 1:atom; j = (atom+1):(2*atom); k=(2*atom+1):ca;
   l = (ca+1):(ca+atom); m=(ca+atom+1):(ca+2*atom); n=(ca+2*atom+1):(cb);

P1=conventional_3x3x6(A(i,i),B(i,i),nmin);
P2=conventional_3x3x6(A(i,i),B(i,j),nmin);
P3=conventional_3x3x6(A(i,i),B(i,k),nmin);
P4=conventional_3x3x6(A(i,i),B(i,l),nmin);
P5=conventional_3x3x6(A(i,i),B(i,m),nmin);
P6=conventional_3x3x6(A(i,i),B(i,n),nmin);
P7=conventional_3x3x6(A(i,j),B(j,i),nmin);
P8=conventional_3x3x6(A(i,j),B(j,j),nmin);
P9=conventional_3x3x6(A(i,j),B(j,k),nmin);
P10=conventional_3x3x6(A(i,j),B(j,l),nmin);
P11=conventional_3x3x6(A(i,j),B(j,m),nmin);
P12=conventional_3x3x6(A(i,j),B(j,n),nmin);
P13=conventional_3x3x6(A(i,k),B(k,i),nmin);
P14=conventional_3x3x6(A(i,k),B(k,j),nmin);
P15=conventional_3x3x6(A(i,k),B(k,k),nmin);
P16=conventional_3x3x6(A(i,k),B(k,l),nmin);
P17=conventional_3x3x6(A(i,k),B(k,m),nmin);
P18=conventional_3x3x6(A(i,k),B(k,n),nmin);
P19=conventional_3x3x6(A(j,i),B(i,i),nmin);
P20=conventional_3x3x6(A(j,i),B(i,j),nmin);
P21=conventional_3x3x6(A(j,i),B(i,k),nmin);
P22=conventional_3x3x6(A(j,i),B(i,l),nmin);
P23=conventional_3x3x6(A(j,i),B(i,m),nmin);
P24=conventional_3x3x6(A(j,i),B(i,n),nmin);
P25=conventional_3x3x6(A(j,j),B(j,i),nmin);
P26=conventional_3x3x6(A(j,j),B(j,j),nmin);
P27=conventional_3x3x6(A(j,j),B(j,k),nmin);
P28=conventional_3x3x6(A(j,j),B(j,l),nmin);
P29=conventional_3x3x6(A(j,j),B(j,m),nmin);
P30=conventional_3x3x6(A(j,j),B(j,n),nmin);
P31=conventional_3x3x6(A(j,k),B(k,i),nmin);
P32=conventional_3x3x6(A(j,k),B(k,j),nmin);
P33=conventional_3x3x6(A(j,k),B(k,k),nmin);
P34=conventional_3x3x6(A(j,k),B(k,l),nmin);
P35=conventional_3x3x6(A(j,k),B(k,m),nmin);
P36=conventional_3x3x6(A(j,k),B(k,n),nmin);
P37=conventional_3x3x6(A(k,i),B(i,i),nmin);
P38=conventional_3x3x6(A(k,i),B(i,j),nmin);
P39=conventional_3x3x6(A(k,i),B(i,k),nmin);
P40=conventional_3x3x6(A(k,i),B(i,l),nmin);
P41=conventional_3x3x6(A(k,i),B(i,m),nmin);
P42=conventional_3x3x6(A(k,i),B(i,n),nmin);
P43=conventional_3x3x6(A(k,j),B(j,i),nmin);
P44=conventional_3x3x6(A(k,j),B(j,j),nmin);
P45=conventional_3x3x6(A(k,j),B(j,k),nmin);
P46=conventional_3x3x6(A(k,j),B(j,l),nmin);
P47=conventional_3x3x6(A(k,j),B(j,m),nmin);
P48=conventional_3x3x6(A(k,j),B(j,n),nmin);
P49=conventional_3x3x6(A(k,k),B(k,i),nmin);
P50=conventional_3x3x6(A(k,k),B(k,j),nmin);
P51=conventional_3x3x6(A(k,k),B(k,k),nmin);
P52=conventional_3x3x6(A(k,k),B(k,l),nmin);
P53=conventional_3x3x6(A(k,k),B(k,m),nmin);
P54=conventional_3x3x6(A(k,k),B(k,n),nmin);

   C = [ 
P1+P7+P13    P2+P8+P14   P3+P9+P15   P4+P10+P16  P5+P11+P17  P6+P12+P18 ;  
P19+P25+P31  P20+P26+P32 P21+P27+P33 P22+P28+P34 P23+P29+P35 P24+P30+P36 ;
P37+P43+P49  P38+P44+P50 P39+P45+P51 P40+P46+P52 P41+P47+P53 P42+P48+P54
	];

end
