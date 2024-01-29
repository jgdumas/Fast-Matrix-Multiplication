function C = DPS_smallrat(A, B, nmin)

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

u1=A(i,j)+A(j,j)*2;
PA1=u1*5/9;
PA2=A(i,i)*8/9-A(i,j)*2/3;
PA3=A(i,i)*4/9+A(j,i)*8/9+u1*2/9;
PA4=A(i,j)*10/9;
PA0=PA2-PA3;
PA5=PA1+PA0;
PA6=PA4+PA0;

v1=B(i,j)/2;
PB1=v1-B(j,j);
PB2=v1-B(i,i);
PB3=B(i,j)*5/4;
PB4=B(j,j)*2/5-B(j,i)*4/5+PB2*3/5;
PB0=PB1+PB4;
PB5=PB0-PB2;
PB6=PB3-PB0;

	d0 = DPS_smallrat(PA0,PB0,nmin);
	d1 = DPS_smallrat(PA1,PB1,nmin);
	d2 = DPS_smallrat(PA2,PB2,nmin);
	d3 = DPS_smallrat(PA3,PB3,nmin);
	d4 = DPS_smallrat(PA4,PB4,nmin);
	d5 = DPS_smallrat(PA5,PB5,nmin);
	d6 = DPS_smallrat(PA6,PB6,nmin);

w1=d6+d0+d4;
w2=d5+d6;
w3=d3+w1;
w4=d2+d4;
w5=d1+w1;
w6=w3*9/20;
c11=w6-w4*9/8;
c12=w3*9/10;
c21=w5*27/40-w2*9/8+c11/2;
c22=w6-w5*9/10;

	C = [ c11 c12 ; c21 c22 ] ;

end ;
