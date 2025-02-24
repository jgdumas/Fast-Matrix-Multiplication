function C = DPS_intermediate_2_2_2(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*2^l x nmin*2^l x nmin*2^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<2)||(k<2)||(n<2)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,2);ku=k-rem(k,2);nu=n-rem(n,2);
  l=ceil(min([log(mu)/log(2),log(ku)/log(2),log(nu)/log(2)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(2)),floor(log(k/nmin)/log(2)),floor(log(n/nmin)/log(2))]);
    mu=nmin*2^l; ku=nmin*2^l; nu=nmin*2^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=DPS_intermediate_2_2_2(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS_intermediate(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS_intermediate(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS_intermediate(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<2;2;2>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
oA3 = A(r0,c1)*289/512-A(r0,c0);
oA2 = A(r1,c1)*334084/345665-A(r0,c1)*51622047/88490240;
oA4 = A(r0,c1)*289/256;
oA6 = (A(r1,c1)+oA3)*167042/345665-A(r1,c0)*295936/345665;
oA0 = oA6-oA2;
oA1 = oA3-oA0;
oA5 = oA0-oA4;

[m,n] = size(B);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
t4 = B(r1,c0)*512/289-B(r0,c0);
oB1 = -B(r0,c0)*345665/295936;
oB2 = t4*345665/591872;
oB3 = -B(r0,c0)*178623/295936-B(r0,c1);
oB6 = B(r0,c1)/2-B(r1,c1)*256/289+t4*289/1024;
oB0 = oB3+oB6;
oB4 = oB2-oB0;
oB5 = oB1-oB0;

iC0 = DPS_intermediate( oA0, oB0, nmin, peeling, level);
iC1 = DPS_intermediate( oA1, oB1, nmin, peeling, level);
iC2 = DPS_intermediate( oA2, oB2, nmin, peeling, level);
iC3 = DPS_intermediate( oA3, oB3, nmin, peeling, level);
iC4 = DPS_intermediate( oA4, oB4, nmin, peeling, level);
iC5 = DPS_intermediate( oA5, oB5, nmin, peeling, level);
iC6 = DPS_intermediate( oA6, oB6, nmin, peeling, level);

b2 = iC0+iC1+iC5;
z4 = iC4+b2;
z3 = iC3-b2;
oC0 = z4*295936/345665;
oC1 = z4*167042/345665+z3;
oC2 = z4*178623/345665+iC2-iC1;
oC3 = (oC2-z3)*289/512-(iC6+iC5)*345665/295936;

C = [ oC0 oC1 ; oC2 oC3 ] ;
  end
end
end
