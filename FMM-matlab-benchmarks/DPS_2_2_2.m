function C = DPS_2_2_2(A, B, nmin, peeling, level)
NTH2ROOT3once=nthroot(3,2);
NTH2ROOT3o2=nthroot(3,2)/2;
NTH2ROOT3t2o=2/nthroot(3,2);
NTH2ROOT3t1o=1/nthroot(3,2);
NTH2ROOT3o4=nthroot(3,2)/4;
NTH2ROOT3o8=nthroot(3,2)/8;
NTH2ROOT3o16=nthroot(3,2)/16;
NTH2ROOT3o3=nthroot(3,2)/3;
NTH2ROOT3f2o3=nthroot(3,2)*2/3;

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
    C(1:mu,1:nu)=DPS_2_2_2(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<2;2;2>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
t4 = A(r1,c0)-A(r0,c0)*NTH2ROOT3once;
r5 = A(r1,c1)*NTH2ROOT3o3;
oA2 = A(r0,c1)+r5;
r7 = A(r1,c1)*NTH2ROOT3once/6;
r8 = (t4-A(r0,c1))/2;
oA0 = A(r0,c0)*NTH2ROOT3o2+(A(r1,c0)+oA2)/2;
oA1 = A(r1,c0)-r5;
oA3 = -A(r1,c1)*NTH2ROOT3f2o3;
oA4 = r8-A(r1,c1)*NTH2ROOT3o2;
oA5 = r7+r8;
oA6 = (A(r0,c1)+t4)/2-r7;

[m,n] = size(B);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
r4 = B(r0,c1)*NTH2ROOT3o3;
oB1 = r4-B(r0,c0);
oB0 = B(r0,c1)*NTH2ROOT3f2o3;
oB2 = r4-B(r1,c1);
oB3 = B(r1,c0)*NTH2ROOT3o2-(oB1+B(r1,c1))/2;
v11 = B(r0,c0)+B(r1,c0)*NTH2ROOT3once;
oB4 = v11-oB3*3-B(r1,c1)*2;
oB6 = B(r0,c0)-B(r1,c1)-oB3;
oB5 = v11-oB3;

iC0 = DPS( oA0, oB0, nmin, peeling, level);
iC1 = DPS( oA1, oB1, nmin, peeling, level);
iC2 = DPS( oA2, oB2, nmin, peeling, level);
iC3 = DPS( oA3, oB3, nmin, peeling, level);
iC4 = DPS( oA4, oB4, nmin, peeling, level);
iC5 = DPS( oA5, oB5, nmin, peeling, level);
iC6 = DPS( oA6, oB6, nmin, peeling, level);
t10 = iC3/2;
t9 = iC0/2;
d4 = iC2+t10;
t5 = (iC5+iC4)/2;
oC0 = iC4*NTH2ROOT3o2-iC6*NTH2ROOT3f2o3+(iC0-iC5)*NTH2ROOT3once/6+(d4-iC1)*NTH2ROOT3o3;
t4 = t9+t5;
oC1 = t9-d4-t5;
oC2 = t4-t10-iC1;
oC3 = iC3*NTH2ROOT3o2+t4*NTH2ROOT3once;

C = [ oC0 oC1 ; oC2 oC3 ] ;
  end
end
end
