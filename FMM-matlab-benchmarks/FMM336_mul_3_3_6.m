function C = FMM336_mul_3_3_6(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*3^l x nmin*3^l x nmin*6^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<3)||(k<3)||(n<6)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,3);ku=k-rem(k,3);nu=n-rem(n,6);
  l=ceil(min([log(mu)/log(3),log(ku)/log(3),log(nu)/log(6)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(3)),floor(log(k/nmin)/log(3)),floor(log(n/nmin)/log(6))]);
    mu=nmin*3^l; ku=nmin*3^l; nu=nmin*6^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMM336_mul_3_3_6(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM336_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM336_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM336_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;3;6>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = A(r2,c2)-A(r0,c2);
t10 = A(r1,c1)-A(r1,c2);
oA22 = A(r0,c0)+A(r1,c2);
oA4 = A(r2,c1)+t10;
oA8 = A(r1,c0)-t9;
t14 = A(r1,c0)-A(r2,c2);
oA24 = -A(r2,c0)-A(r2,c1);
oA27 = A(r1,c1)-A(r2,c0);
oA18 = A(r0,c1)+t9;
oA0 = t14-A(r1,c1);
oA6 = A(r2,c2)+t10;
t20 = A(r2,c1)+A(r1,c1);
oA30 = -A(r0,c1)-A(r2,c2);
t22 = A(r1,c0)+A(r0,c1);
oA1 = A(r1,c0)-t20;
oA5 = A(r2,c1)-t9;
oA9 = oA4-t14;
oA11 = A(r0,c2)+t20;
oA12 = oA8-t10;
oA15 = A(r1,c2)-t9;
oA16 = t22-oA4;
oA19 = oA30-A(r1,c1);
oA20 = A(r0,c0)-oA4;
oA21 = t9-oA22;
oA23 = A(r0,c0)-oA0;
oA25 = oA8-A(r2,c0);
oA26 = oA6-A(r2,c0);
oA28 = -t22;
oA29 = oA18-A(r2,c1);
oA31 = -A(r0,c1)-A(r1,c2);
oA32 = A(r0,c2)+oA22;
oA33 = A(r1,c0)-oA22;
oA35 = A(r2,c1)-oA22;
oA37 = t10-A(r2,c0);
oA38 = -A(r0,c2)-oA27;
oA39 = -oA24-A(r1,c0);
oA2 = A(r2,c2);
oA3 = A(r2,c1);
oA7 = A(r0,c2);
oA10 = A(r1,c2);
oA13 = A(r1,c0);
oA14 = A(r1,c1);
oA17 = A(r0,c1);
oA34 = A(r0,c0);
oA36 = A(r2,c0);

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = B(r1,c0)-B(r0,c4);
t19 = B(r0,c3)-B(r0,c2);
oB13 = B(r1,c2)-t19;
oB6 = B(r0,c3)+t18;
t22 = B(r1,c4)+B(r1,c5);
t23 = B(r2,c1)-B(r2,c2);
t25 = -B(r1,c1)-B(r2,c0);
t26 = B(r0,c1)+oB13;
oB16 = B(r2,c4)+B(r0,c0)-B(r2,c3);
t28 = B(r1,c0)+B(r2,c5);
oB2 = B(r2,c3)-t23;
oB7 = B(r1,c1)+t22;
oB10 = B(r2,c4)+t26;
oB17 = B(r2,c1)+oB13-B(r0,c0);
oB18 = t25-B(r1,c4)-t18;
oB21 = B(r0,c3)+t23+t28;
oB22 = t26-B(r1,c3)-t22;
oB24 = B(r0,c1)-B(r1,c3)-B(r1,c5);
oB25 = B(r2,c2)-t28;
oB26 = B(r2,c3)+B(r2,c5)+oB6;
oB27 = B(r1,c1)-B(r1,c3)-B(r2,c4);
oB29 = B(r0,c5)-B(r1,c4)+B(r2,c1);
oB30 = -B(r0,c5)-B(r1,c1)-B(r2,c3);
oB31 = B(r2,c4)-t18-t19;
oB32 = t22-t23-B(r0,c5);
oB34 = B(r0,c1)+B(r1,c0)+B(r1,c2);
oB36 = t25-oB6;
oB37 = B(r2,c2)-oB16-B(r0,c1);
oB38 = B(r1,c5)-B(r2,c0)-B(r1,c0);
oB0 = B(r1,c0);
oB1 = B(r2,c4);
oB3 = B(r1,c4);
oB4 = B(r0,c1);
oB5 = B(r2,c1);
oB8 = B(r0,c3);
oB9 = B(r2,c3);
oB11 = B(r1,c5);
oB12 = B(r2,c2);
oB14 = B(r1,c1);
oB15 = B(r0,c4);
oB19 = B(r2,c0);
oB20 = B(r1,c3);
oB23 = B(r2,c5);
oB28 = B(r1,c2);
oB33 = B(r0,c2);
oB35 = B(r0,c5);
oB39 = B(r0,c0);

iC0 = FMM336_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMM336_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMM336_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMM336_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMM336_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMM336_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMM336_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMM336_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMM336_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMM336_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMM336_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMM336_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMM336_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMM336_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMM336_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMM336_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMM336_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMM336_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMM336_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMM336_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMM336_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMM336_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMM336_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMM336_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMM336_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMM336_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMM336_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMM336_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMM336_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMM336_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMM336_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMM336_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMM336_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMM336_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMM336_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMM336_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMM336_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMM336_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMM336_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMM336_mul( oA39, oB39, nmin, peeling, level);
b2 = iC37+iC34;
oC17 = iC39+iC36+iC35+iC33;
oC6 = iC31-iC30-iC27+iC26;
t24 = iC29+iC24;
t26 = iC22-iC18;
t28 = iC20+iC16;
b14 = iC29-iC16+iC12;
b17 = iC37-iC27-iC10;
b22 = iC33-iC31-iC23-iC6;
t27 = iC17-iC5;
oC16 = iC23+iC19+t28;
b29 = iC39+t27;
oC9 = iC22-iC35-iC29+iC3+t27;
oC2 = iC21+iC17+t26;
oC7 = iC38+iC32+b2;
oC15 = iC25-iC23-iC12+iC0-b2;
oC11 = iC28+iC25+t24;
oC5 = iC34+iC31+iC15+iC10-t26;
oC14 = iC37+iC24+iC4+b14;
oC0 = iC11-iC32-iC20+b14;
oC12 = iC13-t24+b29;
oC13 = iC8-iC33-iC25-t26-b29;
oC8 = iC7-iC32-iC22-b17;
oC4 = iC2-iC30-iC17+b17;
b31 = iC39+iC27+b22;
oC3 = iC1-iC16+b22;
oC10 = iC9-iC30-b31;
oC1 = iC14-iC35+t28-b31;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
