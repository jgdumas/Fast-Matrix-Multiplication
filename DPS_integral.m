function C = DPS_integral(A, B, nmin)

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

PA1=((13^2)*(A(i,i)+A(j,j))+(7*14)*(A(i,j)-A(j,i)))/(5*17*449);
PA5=A(j,j)/(13^2);
PA3=-((14^2)*A(i,j)+(3*71*89)*PA5)/(5*17*449);
PA4=(PA5-A(j,i)/(2*7^2))/2;
PA2=PA4-PA1;
PA6=PA1-PA5;
PA7=PA1+PA3;
PB1=(7*14)*(B(j,i)-B(i,j))-(13^2)*(B(i,i)+B(j,j));
MB1=-B(i,i)/(13^2);
PB2=(5*17*449)*MB1;
PB3=(5*17*449)*(B(j,i)/(7*14)+MB1)/2;
PB4=(3*71*89)*MB1-(14^2)*B(i,j);
PB5=PB3-PB1;
PB6=PB2-PB1;
PB7=PB1-PB4;
d1=DPS_integral(PA1,PB1,nmin);
d2=DPS_integral(PA2,PB2,nmin);
d3=DPS_integral(PA3,PB3,nmin);
d4=DPS_integral(PA4,PB4,nmin);
d5=DPS_integral(PA5,PB5,nmin);
d6=DPS_integral(PA6,PB6,nmin);
d7=DPS_integral(PA7,PB7,nmin);
PC1=d1+d2+d6+d5;
PC3=PC1/(5*17*449);
PC4=-(3*71*89)*PC3;
PC2=d4+d5;
PC1=d3-d2;
c11=PC4-PC1;
c22=PC2+PC4;
PC1=((2*7)^2)*PC3;
c12=((13)^2)*PC1;
PC2=(5*17*449)*(d6+d7)/(14^2);
c21=((2*7^2)*((c11+c22)-(2*7^2)*PC1)+PC2)/(13^2);


	C = [ c11 c21 ; c12 c22 ] ;

end ;
