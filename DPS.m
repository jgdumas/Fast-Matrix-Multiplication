function C = DPS(A, B, nmin)
%          For A and B are matrices of dimension power of 2,
%                              computes the product C = A*B.
%          Used recursively until dimension <= NMIN, is reached,
%          at which point standard multiplication is used.
%          Growth-factor: 2\sqrt(2)+16\sqrt(3)\approx 12.06603145.
%          Frobenius norms: \sqrt(10),\sqrt(10),\sqrt(10),
%          7+7+10 additions, 4+4+4 multiplications/divisions.

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

SQRT3o2 = sqrt(3)/2;
SQRT3o3 = sqrt(3)/3;

t1 = A(j,j)*SQRT3o3;
PA1 = A(j,i)-t1;
PA2 = A(i,j)+t1;
PA3 = t1*2;
PA0 = A(i,i)*SQRT3o2+(A(j,i)+PA2)/2;
PA4 = PA1-PA0;
PA5 = PA4+PA3;
PA6 = PA4+PA2;

s1 = B(i,j)*SQRT3o3;
PB0 = s1*2;
PB1 = s1-B(i,i);
PB2 = s1-B(j,j);
PB3 = (PB1+B(j,j))/2-B(j,i)*SQRT3o2;
PB4 = PB2+PB3;
PB5 = PB0-PB4;
PB6 = PB4-PB1;

d1 = DPS(PA0,PB0,nmin);
d2 = DPS(PA1,PB1,nmin);
d3 = DPS(PA2,PB2,nmin);
d4 = DPS(PA3,PB3,nmin);
d5 = DPS(PA4,PB4,nmin);
d6 = DPS(PA5,PB5,nmin);
d7 = DPS(PA6,PB6,nmin);

w1 = d7+d6;
w2 = d5+d1+d6;
w3 = w2-d2;
w4 = d4+w2;
w5 = w4/2;
c12 = d1-d3-w5;
c11 = (w3-c12-w1*2)*SQRT3o3;
c21 = w3-w5;
c22 = w4*SQRT3o2;

C = [ c11 c12 ; c21 c22 ] ;
end;
