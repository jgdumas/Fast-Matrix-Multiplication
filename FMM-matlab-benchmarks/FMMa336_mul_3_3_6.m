function C = FMMa336_mul_3_3_6(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMMa336_mul_3_3_6(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMMa336_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMMa336_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMMa336_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;3;6>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = A(r0,c1)-A(r1,c0);
t10 = A(r2,c1)+A(r0,c2);
oA26 = A(r2,c2)-t9;
t12 = A(r0,c0)+A(r0,c1);
t13 = A(r1,c1)+A(r2,c1);
oA7 = A(r1,c1)+t10;
oA12 = A(r0,c0)+t9;
oA4 = A(r0,c1)-t10;
oA31 = A(r2,c1)-A(r2,c0);
oA23 = A(r1,c0)-A(r1,c2);
oA32 = t12-A(r1,c2);
t20 = A(r1,c1)+A(r2,c0);
oA21 = A(r1,c2)-A(r0,c1);
oA0 = t9-A(r0,c2);
oA3 = oA12-t10;
oA5 = A(r1,c1)+t9;
oA10 = A(r0,c0)-t13;
oA11 = t12-A(r2,c1);
oA13 = t12-oA7;
oA14 = A(r1,c0)-t13;
oA17 = t20-A(r0,c0);
oA18 = A(r2,c0)-t10;
oA19 = oA31-A(r1,c0);
oA20 = A(r2,c1)-A(r1,c2);
oA22 = oA32-oA7;
oA24 = A(r1,c0)+A(r2,c2)+t10-A(r0,c0)-A(r0,c1)*2;
oA27 = oA26-t13;
oA28 = -A(r2,c0)-oA4;
oA29 = A(r2,c0)-oA12;
oA30 = -t20;
oA34 = -oA21-A(r0,c2);
oA35 = A(r1,c1)-oA23;
oA36 = A(r0,c1)-A(r2,c2);
oA37 = oA26-A(r0,c0);
oA38 = -oA26-A(r0,c2);
oA39 = A(r1,c1)-oA26;
oA1 = A(r2,c1);
oA2 = A(r0,c0);
oA6 = A(r1,c0);
oA8 = A(r0,c1);
oA9 = A(r1,c1);
oA15 = A(r0,c2);
oA16 = A(r2,c0);
oA25 = A(r2,c2);
oA33 = A(r1,c2);

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = B(r1,c2)-B(r1,c5);
t20 = B(r0,c4)+B(r2,c5);
t21 = B(r1,c1)+B(r1,c3)+t18;
t22 = B(r1,c4)-B(r0,c0);
t23 = B(r2,c4)-B(r0,c5);
t24 = B(r2,c0)+B(r2,c1);
t25 = B(r1,c4)-B(r2,c2);
t26 = B(r0,c1)-B(r0,c4);
oB9 = B(r0,c3)-t21;
t28 = B(r0,c2)-t23;
oB15 = t20-B(r0,c5);
t30 = B(r1,c0)-B(r2,c3);
oB10 = B(r1,c4)+t24;
oB12 = B(r0,c2)+t21;
oB14 = t18-B(r1,c0);
oB16 = B(r2,c1)+t21+t22;
oB19 = -t18-t26-B(r0,c5);
oB20 = B(r2,c0)-B(r2,c3)-B(r1,c5);
oB21 = t20-t28;
oB22 = B(r1,c4)-t30;
oB23 = t23-oB9;
oB25 = B(r0,c3)-B(r2,c5)+t28;
oB27 = B(r1,c2)-t24-t30;
oB28 = B(r0,c4)+t25;
oB31 = B(r2,c1)+B(r2,c2)-B(r0,c5);
oB32 = B(r0,c2)+B(r1,c1)+B(r1,c2);
oB34 = B(r2,c0)+t20+t25;
oB35 = B(r1,c0)+B(r1,c3)-B(r0,c3);
oB36 = B(r1,c0)+t26;
oB37 = B(r0,c2)-t22-t24;
oB38 = B(r0,c1)+B(r1,c2)-oB15;
oB39 = B(r0,c3)+t22;
oB0 = B(r2,c5);
oB1 = B(r2,c1);
oB2 = B(r0,c2);
oB3 = B(r1,c0);
oB4 = B(r2,c0);
oB5 = B(r0,c3);
oB6 = B(r0,c5);
oB7 = B(r1,c2);
oB8 = B(r0,c4);
oB11 = B(r1,c5);
oB13 = B(r1,c4);
oB17 = B(r0,c0);
oB18 = B(r0,c1);
oB24 = B(r2,c3);
oB26 = B(r2,c4);
oB29 = B(r1,c3);
oB30 = B(r1,c1);
oB33 = B(r2,c2);

iC0 = FMMa336_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMMa336_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMMa336_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMMa336_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMMa336_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMMa336_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMMa336_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMMa336_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMMa336_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMMa336_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMMa336_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMMa336_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMMa336_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMMa336_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMMa336_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMMa336_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMMa336_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMMa336_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMMa336_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMMa336_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMMa336_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMMa336_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMMa336_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMMa336_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMMa336_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMMa336_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMMa336_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMMa336_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMMa336_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMMa336_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMMa336_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMMa336_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMMa336_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMMa336_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMMa336_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMMa336_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMMa336_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMMa336_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMMa336_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMMa336_mul( oA39, oB39, nmin, peeling, level);
b1 = iC39+iC33;
t28 = iC37+iC32;
oC9 = iC30-iC31+iC27-iC26;
oC6 = iC29+iC28+iC25+iC24;
oC13 = iC22+iC21-iC18+iC17;
oC2 = iC23+iC20+iC19+iC16;
b16 = iC15-iC38;
b22 = iC8-iC25;
b23 = iC19-iC31-iC6;
b28 = iC38+iC25+iC0;
b30 = iC31+iC18+b16;
oC8 = iC21-iC32-iC26+iC2+b16;
oC7 = iC23-iC26+iC9+iC6-b1;
oC0 = iC38+iC34+t28;
b32 = iC18-iC29+b22;
oC12 = iC36+iC35+b1;
oC3 = iC21+iC5+b22-b1;
oC17 = iC36-iC27+iC14-b23;
oC1 = iC33+iC20+iC1+b23;
b34 = iC29+iC19+b28;
oC14 = iC23-iC32+iC12-b28;
oC15 = iC36+iC3+b32;
oC5 = iC13-iC33-iC24-iC22+b32;
oC11 = iC27+iC7-b30;
oC10 = iC10-iC22-t28+b30;
oC16 = iC11+b34;
oC4 = iC24+iC20+iC4+t28+b34;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
