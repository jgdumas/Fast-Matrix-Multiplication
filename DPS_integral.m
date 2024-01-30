function C = DPS_integral(A, B, nmin)
%          For A and B are matrices of dimension power of 2, computes the product C = A*B.
%          Used recursively until dimension <= NMIN, is reached,
%          at which point standard multiplication is used.
%          Growth-factor:4369701675/473027282+2\sqrt(2) \approx 12.06616423.
%          Frobenius norms:
%                   \frac{3 \sqrt{46527883412970}}{6449885},
%                   \frac{3 \sqrt{46527883412970}}{6492304},
%                   \frac{3 \sqrt{1219124418}}{33124},
%          7+7+11 additions, 7+7+4 multiplications/divisions.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

t1=A(i,j)+A(j,j)*18957/33124;
PA1=A(i,i)*33124/38165+A(j,i)*18957/38165+t1*19208/38165;
PA2=t1*38416/38165;
PA3=A(j,i)-A(j,j)*98/169;
PA4=A(j,j)*196/169;
PA0=PA1-PA3;
PA5=PA0-PA4;
PA6=PA0-PA2;

s1=B(i,i)-B(j,i)*169/98;
PB1=B(i,i)*38165/33124;
PB2=s1*38165/66248;
PB3=B(i,i)*18957/33124+B(i,j);
PB6=B(i,j)/2-B(j,j)*169/196-s1*49/169;
PB0=PB6-PB3;
PB4=PB2+PB0;
PB5=PB1+PB0;

d0=DPS_integral(PA0,PB0,nmin);
d1=DPS_integral(PA1,PB1,nmin);
d2=DPS_integral(PA2,PB2,nmin);
d3=DPS_integral(PA3,PB3,nmin);
d4=DPS_integral(PA4,PB4,nmin);
d5=DPS_integral(PA5,PB5,nmin);
d6=DPS_integral(PA6,PB6,nmin);

w0=d1-d4-d5;
w1=d6-d5;
w2=d1-d2;
w3=d3-d4;
w4=w3+w2;
w5=d0+w0;
w6=w5*18957/38165;
c11=w2-w6;
c21=w5*33124/38165;
c12=w1*38165/33124-c21+w4*98/169;
c22=w3-w6;


	C = [ c11 c12 ; c21 c22 ] ;

end ;
