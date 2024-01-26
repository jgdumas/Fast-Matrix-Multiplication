function C = DPS_evenpow(A, B, nmin)

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
error('The matrix dimension must be a power of 2.')
end

if n <= nmin
C = A*B;
else
m = n/2; i = 1:m; j = m+1:n;

CSTE=A(i,i)/2;
PA0=-5*CSTE;
PA6=-5*(CSTE+A(j,i))/2;
PA4=3*CSTE-2*A(i,j);
PA1=2*(A(i,i)+A(j,j))-A(i,j)+A(j,i);
PA2=PA1+PA6;
PA3=PA1-PA4;
PA5=PA1+PA0;
CSTE2=B(j,i)/2;
PB0=5*CSTE2/2;
PB1=(4*(B(i,j)+B(j,i))+3*(B(j,j)-B(i,i)))/5;
PB4=CSTE2+B(j,j);
PB6=-B(i,i)+CSTE2;
PB2=PB4-PB1;
PB5=PB0-PB1;
PB7=PB1-PB6;
d1=DPS_evenpow(PA1,PB1,nmin);
d2=DPS_evenpow(PA2,PB2,nmin);
d3=DPS_evenpow(PA3,PB0,nmin);
d4=DPS_evenpow(PA4,PB4,nmin);
d5=DPS_evenpow(PA5,PB5,nmin);
d6=DPS_evenpow(PA6,PB6,nmin);
d7=DPS_evenpow(PA0,PB7,nmin);
PC1=(d2-d5)/2;
PC2=-(d4+d7)/2;
PC3=-2*(d3+d6)/5+PC2;
c11=2*(-d1+d3-d5-d7)/5;
PC4=c11/2;
c21=PC2-PC4;
c12=c21-PC3;
c22=PC1+(c21+3*(c12-PC4)/2)/2;

C = [ c11 c21 ; c12 c22 ] ;

end ;
