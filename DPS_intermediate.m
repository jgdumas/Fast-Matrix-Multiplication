function C = DPS_intermediate(A, B, nmin)

if nargin < 3, nmin = 8; end

n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end

if n <= nmin
   C = A*B;
else
	m = n/2; i = 1:m; j = m+1:n;

PA1=(295936*(A(i,j)-A(j,i))-167042*(A(i,i)+A(j,j)))/345665;
PA5=A(i,j)*289/256;
PA3=A(j,j)*334084/345665-A(i,j)*51622047/88490240;
PA4=PA5/2-A(i,i);
PA2=PA4-PA1;
PA6=PA1-PA5;
PA7=PA1+PA3;
PB1=(B(j,i)-B(i,j))/2-(B(i,i)+B(j,j))*256/289;
PB2=-B(i,i)*345665/295936;
PB3=PB2/2+B(j,i)*345665/334084;
PB4=-B(i,i)*178623/295936-B(i,j);
PB5=PB3-PB1;
PB6=PB2-PB1;
PB7=PB1-PB4;
d1=DPS_intermediate(PA1,PB1,nmin);
d2=DPS_intermediate(PA2,PB2,nmin);
d3=DPS_intermediate(PA3,PB3,nmin);
d4=DPS_intermediate(PA4,PB4,nmin);
d5=DPS_intermediate(PA5,PB5,nmin);
d6=DPS_intermediate(PA6,PB6,nmin);
d7=DPS_intermediate(PA7,PB7,nmin);
PC1=d4+d5;
PC2=d3-d2;
PC3=d1+d2+d5+d6;
c11=PC3*295936/345665;
MC1=c11*178623/295936;
c21=PC1-MC1;
c12=PC2+MC1;
c22=(289*(289*c11/512+(c12-c21))-345665*(d6+d7)/(2*289))/512;

C = [ c11 c21 ; c12 c22 ] ;

end ;
