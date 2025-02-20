function C = DPS_smallrat_2_2_2(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
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
    C(1:mu,1:nu)=DPS_smallrat_2_2_2(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS_smallrat(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS_smallrat(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS_smallrat(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<2;2;2>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
t4 = A(r0,c1)+A(r1,c1)*2;
oA1 = t4*5/9;
oA2 = A(r0,c0)*8/9-A(r0,c1)*2/3;
oA3 = A(r0,c0)*4/9+A(r1,c0)*8/9+t4*2/9;
oA4 = -A(r0,c1)*10/9;
oA0 = oA2-oA3;
oA5 = oA1+oA0;
oA6 = oA4-oA0;

[m,n] = size(B);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
r4 = B(r0,c1)/2;
oB2 = r4-B(r0,c0);
oB1 = r4-B(r1,c1);
oB3 = B(r0,c1)*5/4;
oB4 = B(r1,c0)*4/5-B(r1,c1)*2/5-oB2*3/5;
oB0 = oB1-oB4;
oB5 = oB0-oB2;
oB6 = oB0-oB3;

iC0 = DPS_smallrat( oA0, oB0, nmin, peeling, level);
iC1 = DPS_smallrat( oA1, oB1, nmin, peeling, level);
iC2 = DPS_smallrat( oA2, oB2, nmin, peeling, level);
iC3 = DPS_smallrat( oA3, oB3, nmin, peeling, level);
iC4 = DPS_smallrat( oA4, oB4, nmin, peeling, level);
iC5 = DPS_smallrat( oA5, oB5, nmin, peeling, level);
iC6 = DPS_smallrat( oA6, oB6, nmin, peeling, level);

b1 = iC0+iC6+iC4;
z3 = iC3+b1;
z1 = iC1+b1;
t5 = z3*9/20;
oC1 = z3*9/10;
oC0 = t5-(iC2+iC4)*9/8;
oC3 = t5-z1*9/10;
oC2 = z1*27/40-(iC5+iC6)*9/8+oC0/2;

C = [ oC0 oC1 ; oC2 oC3 ] ;
  end
end
end
