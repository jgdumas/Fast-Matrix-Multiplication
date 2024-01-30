function C = DPS_intermediate(A, B, nmin)
%          For A and B are matrices of dimension power of 2, computes the product C = A*B.
%          Used recursively until dimension <= NMIN, is reached,
%          at which point standard multiplication is used.
%          Growth-factor:358452876675/38788923392+2\sqrt(2)\approx 12.06954148.
%          Frobenius norms:
%                    \sqrt(302819615879344530)/176980480 \approx 3.109328685
%                    \sqrt(302819615879344530)/171051008 \approx 3.217113361
%                    \sqrt( 876049400082)/295936 \approx 3.162761904
%          7+7+10 additions, 6+7+5 multiplications/divisions.

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

PA2=A(j,j)*334084/345665-A(i,j)*51622047/88490240;
PA3=A(i,i)-A(i,j)*289/512;
t1=A(j,j)-PA3;
PA4=A(i,j)*289/256;
PA6=t1*167042/345665-A(j,i)*295936/345665;
PA0=PA6-PA2;
PA1=PA0+PA3;
PA5=PA4-PA0;

t4=B(j,i)*512/289-B(i,i);
PB1=B(i,i)*345665/295936;
PB2=t4*345665/591872;
PB3=B(i,i)*178623/295936+B(i,j);
PB6=B(i,j)/2-B(j,j)*256/289+t4*289/1024;
PB0=PB6-PB3;
PB5=PB1+PB0;
PB4=PB2-PB0;

d0=DPS_intermediate(PA0,PB0,nmin);
d1=DPS_intermediate(PA1,PB1,nmin);
d2=DPS_intermediate(PA2,PB2,nmin);
d3=DPS_intermediate(PA3,PB3,nmin);
d4=DPS_intermediate(PA4,PB4,nmin);
d5=DPS_intermediate(PA5,PB5,nmin);
d6=DPS_intermediate(PA6,PB6,nmin);

w1=d5+d0+d4;
w2=d6+d5;
w3=d2+w1;
w4=d1+w1;
c11=w4*295936/345665;
c12=d3+d4-w4*178623/345665;
c21=w3-w4*167042/345665;
c22=(w3-c12)*289/512-w2*345665/295936;

C = [ c11 c12 ; c21 c22 ] ;

end ;
