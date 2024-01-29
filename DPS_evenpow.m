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

t1=A(j,j)/2;
PA2=A(j,i)-t1;
t3=PA2-A(i,j);
PA3=A(i,j)-t1;
PA4=A(j,i)+t1;
PA5=A(i,i)+t3/2;
PA0=PA2-PA3;
PA6=PA4-PA0;
PA1=PA5-PA0;

s1=B(i,j)/2;
PB1=B(i,i)+s1;
s3=PB1-B(j,j);
PB2=s1-B(j,j);
PB3=s3/2-B(j,i);
PB4=B(j,j)+s1;
PB0=PB1-PB4;
PB5=PB0-PB2;
PB6=PB0-PB3;

	d0 = DPS_evenpow(PA0,PB0,nmin);
	d1 = DPS_evenpow(PA1,PB1,nmin);
	d2 = DPS_evenpow(PA2,PB2,nmin);
	d3 = DPS_evenpow(PA3,PB3,nmin);
	d4 = DPS_evenpow(PA4,PB4,nmin);
	d5 = DPS_evenpow(PA5,PB5,nmin);
	d6 = DPS_evenpow(PA6,PB6,nmin);

r1=d4+d2;
r2=d4-d2;
r3=d6+d5-d3+d1;
r4=d0+r2/2;
c11=r3/2+r1/4;
c12=d1-d5+r4;
c21=d3+d6+r4;
c22=r1;

	C = [ c11 c12 ; c21 c22 ] ;

end ;
