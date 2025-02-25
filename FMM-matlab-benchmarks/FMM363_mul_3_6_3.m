function C = FMM363_mul_3_6_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*3^l x nmin*6^l x nmin*3^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<3)||(k<6)||(n<3)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,3);ku=k-rem(k,6);nu=n-rem(n,3);
  l=ceil(min([log(mu)/log(3),log(ku)/log(6),log(nu)/log(3)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(3)),floor(log(k/nmin)/log(6)),floor(log(n/nmin)/log(3))]);
    mu=nmin*3^l; ku=nmin*6^l; nu=nmin*3^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMM363_mul_3_6_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM363_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM363_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM363_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;6;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = A(r0,c4)-A(r0,c5);
oA11 = A(r2,c3)-t18;
t20 = A(r1,c4)+A(r2,c4);
t21 = A(r1,c1)+A(r2,c0);
t22 = A(r0,c0)-A(r0,c2);
t23 = A(r1,c5)-A(r2,c5);
t24 = A(r2,c1)+A(r2,c5);
t26 = oA11-A(r0,c3);
t27 = t20+A(r1,c3)+A(r1,c2);
t28 = A(r0,c1)-A(r2,c2);
oA2 = t21-A(r0,c1);
oA8 = A(r0,c5)+t20;
oA9 = A(r1,c3)-t22;
oA12 = oA11+t23;
oA16 = t22-t26-A(r1,c5);
oA17 = A(r0,c1)+A(r1,c0)+A(r1,c4);
oA18 = A(r0,c5)-A(r1,c0)-A(r2,c0);
oA20 = A(r0,c0)-t26;
oA21 = A(r1,c0)+t20+t21;
oA23 = A(r0,c3)+A(r1,c3)-A(r2,c5);
oA25 = t18-t23-A(r2,c4);
oA26 = A(r2,c2)-A(r1,c1)-A(r1,c3);
oA27 = t28-A(r0,c2);
oA28 = A(r0,c4)+A(r1,c4)-A(r1,c5);
oA30 = t22+t28-A(r2,c0);
oA32 = t24-oA11-t21;
oA35 = t27-A(r0,c0);
oA36 = A(r0,c5)+t27;
oA37 = A(r1,c5)+A(r2,c1)-A(r0,c1);
oA38 = t24-A(r1,c1);
oA39 = A(r1,c2)+A(r1,c4)-A(r0,c2);
oA0 = A(r2,c5);
oA1 = A(r0,c2);
oA3 = A(r0,c5);
oA4 = A(r1,c5);
oA5 = A(r2,c4);
oA6 = A(r1,c3);
oA7 = A(r2,c0);
oA10 = A(r0,c1);
oA13 = A(r1,c4);
oA14 = A(r0,c0);
oA15 = A(r1,c1);
oA19 = A(r0,c3);
oA22 = A(r1,c0);
oA24 = A(r0,c4);
oA29 = A(r2,c3);
oA31 = A(r2,c2);
oA33 = A(r1,c2);
oA34 = A(r2,c1);

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
oB5 = B(r1,c1)+B(r3,c2)-B(r0,c2);
t20 = -B(r2,c0)-B(r0,c1);
t21 = B(r5,c1)-B(r4,c2);
t22 = B(r5,c2)+B(r3,c1);
t23 = B(r1,c0)+B(r1,c2);
oB14 = t20-B(r4,c0);
t25 = oB5-B(r5,c0);
oB2 = B(r2,c0)-t25;
oB11 = t21-oB14;
oB10 = B(r0,c2)+t23;
t29 = B(r3,c0)-B(r0,c2);
oB0 = t22-B(r2,c2);
t31 = B(r2,c0)-B(r3,c2);
t32 = B(r4,c1)+B(r5,c2);
oB16 = B(r1,c0)-t31;
oB18 = t32-B(r5,c1);
oB19 = oB11-oB0-B(r4,c1);
oB20 = oB14-B(r1,c0)-B(r2,c1);
oB22 = B(r2,c1)-B(r5,c1)+oB10;
oB23 = B(r0,c0)-t22-t25;
oB24 = B(r2,c1)-t21+t23;
oB25 = B(r2,c2)+oB5-B(r0,c0);
oB26 = B(r0,c0)-B(r5,c2)+oB2;
oB28 = B(r2,c2)-t29;
oB29 = B(r4,c0)-B(r4,c2)+oB5;
oB31 = B(r1,c0)+B(r3,c0)-B(r3,c1);
oB32 = B(r5,c1)-t20-t25;
oB34 = B(r1,c2)+t22-t29;
oB36 = t32-B(r2,c2)-t21;
oB37 = B(r5,c0)-t23+t31;
oB1 = B(r1,c0);
oB3 = B(r4,c2);
oB4 = B(r1,c2);
oB6 = B(r3,c1);
oB7 = B(r5,c1);
oB8 = B(r2,c2);
oB9 = B(r2,c0);
oB12 = B(r5,c0);
oB13 = B(r0,c2);
oB15 = B(r5,c2);
oB17 = B(r1,c1);
oB21 = B(r0,c0);
oB27 = B(r2,c1);
oB30 = B(r0,c1);
oB33 = B(r3,c0);
oB35 = B(r4,c0);
oB38 = B(r4,c1);
oB39 = B(r3,c2);

iC0 = FMM363_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMM363_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMM363_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMM363_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMM363_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMM363_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMM363_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMM363_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMM363_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMM363_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMM363_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMM363_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMM363_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMM363_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMM363_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMM363_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMM363_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMM363_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMM363_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMM363_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMM363_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMM363_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMM363_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMM363_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMM363_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMM363_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMM363_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMM363_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMM363_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMM363_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMM363_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMM363_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMM363_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMM363_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMM363_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMM363_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMM363_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMM363_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMM363_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMM363_mul( oA39, oB39, nmin, peeling, level);
oC6 = iC37-iC39-iC38-iC36+iC27+iC26+iC25+iC24;
t18 = iC32+iC23;
b13 = iC30+iC19;
t19 = iC19-iC14;
t15 = iC33+iC8;
t17 = iC36+iC6;
b31 = iC31-iC27-iC10+iC1;
b33 = iC38+iC0;
b34 = iC25-iC23+iC12+iC0;
b35 = iC11-iC20-t19;
oC4 = iC35-iC34+iC33-iC22+iC21-iC20-t18;
oC3 = iC29-iC31-iC28+iC18+iC17+iC16-b13;
oC1 = iC26+iC11+iC9+iC7+iC6+iC2+t18-b13-b33;
b43 = iC39-iC16+iC9-iC1-b34;
b46 = iC7-iC22-iC18+t15+b31;
oC8 = iC15-iC34-b33-b31+b35;
oC2 = iC13-iC28+t15-b43;
oC0 = iC4-iC24+t17+b35+b43;
oC7 = iC3-t17-b46;
oC5 = iC35-iC29+iC5+t19+b34+b46;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
