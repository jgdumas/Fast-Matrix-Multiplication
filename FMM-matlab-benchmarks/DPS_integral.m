function C = DPS_integral(A, B, nmin)
%          For A and B are matrices of dimension power of 2, computes the product C = A*B.
%          Used recursively until dimension <= NMIN, is reached,
%          at which point standard multiplication is used.
%          Growth-factor:4369701675/473027282+2\sqrt(2) \approx 12.06616423.
%          Frobenius norms:
%                   \frac{3 \sqrt{46527883412970}}{6449885},
%                   \frac{3 \sqrt{46527883412970}}{6492304},
%                   \frac{3 \sqrt{1219124418}}{33124},
%          7+7+10 additions, 6+7+5 multiplications/divisions.

%          Reference:
%          J-G. Dumas, C. Pernet, A. Sedoglavic
%          Strassen's algorithm is not optimally accurate, Feb. 2024
%          https://hal.science/hal-04441653

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

PA2=A(i,j)*38416/38165+A(j,j)*3715572/6449885;
PA3=A(j,i)-A(j,j)*98/169;
PA4=A(j,j)*196/169;
PA6=A(i,i)*33124/38165-(A(i,j)+PA3)*19208/38165;
PA0=PA2+PA6;
PA1=PA3+PA0;
PA5=PA0-PA4;

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

w0=d0+d1-d5;
w1=d1-d2;
w2=w0-d4;
w3=d3-w0;
w4=d6-d5;
c11=w1-w2*18957/38165;
c21=w2*33124/38165;
c22=w2*19208/38165+w3;
w5=w3+c11;
c12=w4*38165/33124+w5*98/169;

	C = [ c11 c12 ; c21 c22 ] ;

end
