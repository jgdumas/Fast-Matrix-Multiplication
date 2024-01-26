function C = DPS(A, B, nmin)

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

r4 = A(j,j)*SQRT3o3;
t4 = A(i,j)+r4;
t5 = A(j,i)+t4;
PA0 = A(i,i)*SQRT3o2+t5/2;
PA1 = A(j,i)-r4;
PA2 = t4;
PA3 = 2*r4;
PA4 = PA1-PA0;
PA5 = PA4+PA3;
PA6 = PA4+PA2;

r4 = B(i,j)*SQRT3o3;
t4 = r4-B(i,i);
t5 = t4+B(j,j);
PB0 = 2*r4;
PB1 = t4;
PB2 = r4-B(j,j);
PB3 = t5/2-B(j,i)*SQRT3o2;
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

t6 = d7+d6;
t4 = d5+d1+d6;
t1 = t4-d2;
t3 = d4+t4;
r3 = t3/2;
c21 = d1-d3-r3;
c12 = t1-r3;
c22 = SQRT3o2*t3;
c11 = SQRT3o3*(t1-c21-2*t6);

C = [ c11 c21 ; c12 c22 ] ;
end;
