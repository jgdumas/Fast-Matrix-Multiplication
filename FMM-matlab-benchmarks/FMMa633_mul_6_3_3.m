function C = FMMa633_mul_6_3_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*6^l x nmin*3^l x nmin*3^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<6)||(k<3)||(n<3)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,6);ku=k-rem(k,3);nu=n-rem(n,3);
  l=ceil(min([log(mu)/log(6),log(ku)/log(3),log(nu)/log(3)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(6)),floor(log(k/nmin)/log(3)),floor(log(n/nmin)/log(3))]);
    mu=nmin*6^l; ku=nmin*3^l; nu=nmin*3^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMMa633_mul_6_3_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMMa633_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMMa633_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMMa633_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<6;3;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t18 = A(r3,c0)+A(r3,c2);
t19 = A(r2,c1)+A(r0,c1);
oA3 = t19-A(r4,c1);
t21 = A(r1,c2)-A(r3,c1);
t22 = A(r0,c0)+A(r2,c0);
oA10 = A(r1,c1)+t18;
oA5 = t22-A(r4,c2);
oA30 = -A(r2,c0)-A(r4,c0)-A(r5,c2);
t27 = A(r3,c0)+A(r0,c2);
t28 = A(r5,c0)-t21;
t29 = A(r1,c0)-oA3;
t30 = A(r0,c0)-A(r5,c1);
t31 = A(r2,c1)-A(r2,c2);
oA6 = A(r2,c1)-t21;
oA11 = t29-A(r4,c0);
oA17 = t18+t30-A(r4,c2);
oA19 = t21-A(r4,c0)-t19;
oA20 = A(r4,c0)-t27;
oA21 = t31+oA5;
oA22 = A(r0,c2)-A(r1,c0)+oA10;
oA23 = A(r0,c0)-A(r2,c2)-A(r3,c1);
oA24 = A(r1,c1)+t27-t29;
oA26 = t22+t31-A(r1,c2);
oA28 = t18-t28;
oA29 = A(r5,c2)-oA3+oA5;
oA32 = A(r1,c0)+A(r4,c2)+oA30;
oA33 = A(r2,c1)-A(r3,c0)+t28;
oA34 = A(r1,c2)-A(r5,c0)+oA10;
oA37 = t30-A(r1,c1);
oA38 = A(r1,c0)-A(r1,c2)+A(r4,c1);
oA39 = A(r2,c0)-A(r3,c0)+A(r5,c1);
oA0 = A(r3,c1);
oA1 = A(r3,c0);
oA2 = A(r4,c2);
oA4 = A(r1,c1);
oA7 = A(r1,c0);
oA8 = A(r2,c1);
oA9 = A(r2,c0);
oA12 = A(r0,c0);
oA13 = A(r3,c2);
oA14 = A(r4,c0);
oA15 = A(r1,c2);
oA16 = A(r5,c1);
oA18 = A(r4,c1);
oA25 = A(r2,c2);
oA27 = A(r0,c2);
oA31 = A(r5,c0);
oA35 = A(r5,c2);
oA36 = A(r0,c1);

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
oB30 = -B(r0,c0)-B(r1,c0);
oB2 = -oB30-B(r1,c0);
v41 = B(r2,c0)-B(r0,c2);
v42 = B(r1,c2)-B(r2,c0);
oB7 = B(r0,c1)+v41;
v44 = -B(r1,c1)-B(r0,c2);
oB9 = B(r1,c1)-v42;
v46 = oB2-oB7;
oB11 = oB2-v44;
v48 = oB30+oB7;
oB38 = -B(r2,c2)-B(r0,c1);
oB29 = B(r1,c0)-B(r1,c1);
oB27 = B(r2,c2)-v41;
oB31 = oB30+v41;
oB10 = oB2-v41;
oB25 = B(r2,c2)+v42;
oB1 = B(r1,c2)+v44;
oB13 = B(r1,c2)+v46;
oB0 = v42-B(r0,c1);
oB28 = v48-B(r1,c2);
oB39 = oB9-B(r2,c2);
oB20 = v44-B(r2,c1);
oB16 = -oB9-oB30;
oB35 = B(r1,c1)+B(r2,c1);
oB32 = oB2-B(r2,c1);
oB36 = B(r2,c0)-B(r2,c2);
oB33 = B(r1,c2)+B(r2,c1);
oB19 = oB30-B(r0,c2);
oB34 = -B(r0,c1)-B(r2,c1);
oB24 = -oB11-oB38;
oB3 = B(r1,c1)+v46;
oB23 = -v42-B(r2,c1);
oB4 = -v44-B(r0,c1);
oB12 = oB2+v42;
oB37 = B(r2,c2)-oB2;
oB22 = v46-B(r2,c1);
oB18 = -v48;
oB5 = B(r1,c1);
oB6 = B(r2,c0);
oB8 = B(r1,c2);
oB14 = B(r0,c2);
oB15 = B(r0,c1);
oB17 = B(r1,c0);
oB21 = B(r2,c1);
oB26 = B(r2,c2);

iC0 = FMMa633_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMMa633_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMMa633_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMMa633_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMMa633_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMMa633_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMMa633_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMMa633_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMMa633_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMMa633_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMMa633_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMMa633_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMMa633_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMMa633_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMMa633_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMMa633_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMMa633_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMMa633_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMMa633_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMMa633_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMMa633_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMMa633_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMMa633_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMMa633_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMMa633_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMMa633_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMMa633_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMMa633_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMMa633_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMMa633_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMMa633_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMMa633_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMMa633_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMMa633_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMMa633_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMMa633_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMMa633_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMMa633_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMMa633_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMMa633_mul( oA39, oB39, nmin, peeling, level);
t34 = iC38+iC34;
t28 = iC36+iC33;
t30 = iC30-iC26;
b5 = iC21-iC38;
t31 = iC23+iC19;
b9 = iC14-iC19;
t33 = iC32-iC12;
oC7 = iC27-iC22+iC10+iC7+t34;
t23 = iC29-t33;
b22 = iC0-iC23+t33;
oC4 = iC15-iC32-iC26+iC2+b5;
oC17 = iC20+iC16+t31;
oC14 = iC37+iC32+t34;
b27 = iC7-iC15-t30;
oC15 = iC31-iC27-t30;
oC0 = iC39+iC35+t28;
oC11 = iC13-iC24-iC22-iC3-t28;
b28 = b9-iC23;
oC12 = iC36+iC6+t30+b9;
oC2 = iC35-iC30+iC9-b28;
oC5 = iC1-iC27-iC16+t28+b28;
oC16 = iC11+t31+t23;
oC10 = iC24-iC16+iC4-t34+t23;
b29 = iC38+b27;
b30 = iC3-iC29+b27;
b31 = iC38+b22;
oC6 = iC18-b29;
oC3 = iC22+iC21+iC17-b29;
oC9 = iC25+b31;
oC13 = iC29+iC28+iC24-b31;
oC8 = iC35+iC5+b5-b30;
oC1 = iC38*2+iC36+iC8+b22+b30;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
