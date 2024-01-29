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

	t1 = (A(j,i)-A(i,j));
	t2 = A(j,j)/2;
	t3 = t1/2;
	t4 = t2/2-A(i,i);
	s1 = -(B(i,j)+B(j,i));
	s2 = -B(i,i)/2;
	s3 = s1/2;
	s4 = s2/2-B(j,j);
	d1 = DPS_evenpow(t1,s1,nmin);
	d2 = DPS_evenpow((t4-t3),(s2-B(i,j)),nmin);
	d3 = DPS_evenpow((A(i,j)-t2),(s4-s3),nmin);
	d4 = DPS_evenpow((t2-A(j,i)),(s2+B(j,i)),nmin);
	d5 = DPS_evenpow((A(i,j)+t2),(s3+s4),nmin);
	d6 = DPS_evenpow((t3+t4),(s2+B(i,j)),nmin);
	d7 = DPS_evenpow((A(j,i)+t2),(B(j,i)-s2),nmin);
	r8 = (d4-d7)/2-d1;
	c11 = d2+d6-r8;
	c12 = d4+d7;
	c21 = (d2-d3-d5-d6)/2-c12/4;
	c22 = d3-d5+r8;

	C = [ c11 c21 ; c12 c22 ] ;

end ;
