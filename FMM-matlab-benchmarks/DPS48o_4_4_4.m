function C = DPS48o_4_4_4(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*4^l x nmin*4^l x nmin*4^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<4)||(k<4)||(n<4)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,4);ku=k-rem(k,4);nu=n-rem(n,4);
  l=ceil(min([log(mu)/log(4),log(ku)/log(4),log(nu)/log(4)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(4)),floor(log(k/nmin)/log(4)),floor(log(n/nmin)/log(4))]);
    mu=nmin*4^l; ku=nmin*4^l; nu=nmin*4^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=DPS48o_4_4_4(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS48o(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS48o(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS48o(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<4;4;4>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
x32 = A(r3,c3)-A(r1,c3);
x37 = A(r3,c3)+A(r1,c3);
x34 = A(r0,c3)+A(r2,c3);
x26 = A(r0,c3)-A(r2,c3);
x31 = A(r1,c1)-A(r3,c1);
x30 = A(r1,c1)+A(r3,c1);
x29 = A(r2,c2)-A(r0,c2);
x28 = A(r2,c2)+A(r0,c2);
x40 = A(r0,c1)+A(r0,c2);
x41 = A(r2,c1)+A(r2,c2);
x36 = A(r3,c1)+A(r3,c2);
x35 = A(r1,c1)+A(r1,c2);
x33 = A(r3,c0)+A(r1,c0);
oA2 = A(r3,c0)-A(r1,c0)-x31;
x27 = A(r2,c0)-A(r0,c0);
oA9 = A(r2,c0)+A(r0,c0)-x28;
x38 = x33-x37;
x39 = x26+x27;
oA15 = x27-x29;
oA11 = x34+x28;
oA7 = x34-x28;
oA34 = oA11-x41;
oA26 = x40-oA34;
oA5 = x30+x33;
oA13 = x30-x33;
oA23 = x31+x32;
oA28 = x31-x32;
oA24 = x26+x29;
oA25 = x30+x37;
oA38 = oA23-x35;
x22 = oA13-oA23;
oA4 = x27+x29;
oA37 = oA4-x41;
x17 = oA15-oA5;
x23 = oA24-oA25;
oA44 = x17-x23;
oA47 = x17+x23;
x20 = oA7-oA28;
x19 = oA2-oA9;
x42 = oA2+oA28;
x43 = oA7-oA9;
x25 = x42+x43;
x24 = x42-x43;
x18 = oA37+oA34;
oA41 = oA13-x35;
oA8 = x36-oA41;
x16 = oA38+oA41;
oA36 = x22+x18;
oA32 = x22-x18;
oA40 = x20+x23;
oA43 = x20-x23;
x21 = oA4-oA11;
oA45 = x21+x16;
oA42 = x21-x16;
oA39 = x17-x19;
oA33 = x17+x19;
oA46 = x19-x20;
oA35 = x19+x20;
oA6 = x36+oA38;
oA10 = x24-oA32;
oA14 = x24-oA42;
oA12 = oA26+x18;
oA1 = oA45+x25;
oA31 = oA36-x25;
x15 = x39-x38;
x10 = x39+x38;
oA21 = x15-oA45;
oA18 = x15-oA36;
oA20 = x10+oA32;
oA27 = x10-oA42;
x14 = x35*2;
x13 = x36*2;
x12 = x41*2;
x11 = x40*2;
oA16 = x14+oA33;
oA17 = x14+oA40;
oA0 = x11+oA40;
oA19 = x11-oA33;
oA3 = oA39+x13;
oA30 = oA43-x13;
oA29 = oA43+x12;
oA22 = oA39+x12;

[m,n] = size(B);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
oB19 = B(r1,c1)-B(r1,c3);
oB12 = B(r1,c3)-B(r1,c0);
oB0 = B(r1,c1)+B(r1,c3);
oB26 = B(r1,c1)+B(r1,c2);
oB24 = B(r3,c3)-B(r3,c0);
oB25 = B(r3,c3)+B(r3,c0);
oB2 = B(r0,c1)-B(r0,c2);
oB9 = B(r0,c1)+B(r0,c2);
oB14 = B(r0,c1)+B(r0,c3);
oB3 = B(r2,c0)+B(r2,c2);
oB30 = B(r2,c0)-B(r2,c2);
oB6 = B(r2,c2)-B(r2,c1);
oB8 = B(r2,c0)+B(r2,c3);
oB10 = B(r0,c0)-B(r0,c2);
oB18 = B(r3,c0)+B(r3,c2);
oB21 = B(r3,c3)-B(r3,c1);
y12 = (B(r1,c2)-B(r0,c2)-B(r2,c2)-B(r3,c2))/2;
oB34 = oB26-y12;
oB36 = oB18+y12;
y14 = (B(r0,c3)+B(r2,c3)+B(r3,c3)-B(r1,c3))/2;
oB41 = y14-oB8;
oB42 = oB14-y14;
y18 = (B(r0,c0)+B(r2,c0)+B(r3,c0)-B(r1,c0))/2;
oB37 = y18-oB12;
oB4 = y14-oB37;
oB32 = oB10-y18;
y21 = oB42-oB32;
oB15 = y21-oB9;
y17 = oB32+oB42;
oB5 = oB2-y17;
oB20 = y12-oB32;
oB13 = oB41+y18;
oB31 = oB36-y18;
y19 = (B(r1,c1)-B(r0,c1)-B(r3,c1)-B(r2,c1))/2;
oB38 = y19-oB6;
oB23 = y12-oB38;
y10 = oB41+oB38;
oB16 = y10+oB3;
oB45 = oB21-y19;
oB1 = oB45-y14;
oB11 = y19-oB34;
oB27 = oB42+y19;
y11 = oB45-oB36;
y15 = oB45+oB36;
y16 = oB34-oB37;
y20 = oB34+oB37;
oB22 = y20-oB19;
oB28 = y15-oB25;
oB29 = oB0-y16;
oB7 = oB24-y11;
y13 = oB41-oB38;
oB17 = y13+oB30;
oB43 = (y13+y16)/2;
oB40 = y16-oB43;
oB33 = (y10+y20)/2;
oB39 = oB33-y20;
oB35 = (y17+y15)/2;
oB47 = oB35-y15;
oB46 = (y11+y21)/2;
oB44 = oB46-y21;

iC0 = DPS48o( oA0, oB0, nmin, peeling, level);
iC1 = DPS48o( oA1, oB1, nmin, peeling, level);
iC2 = DPS48o( oA2, oB2, nmin, peeling, level);
iC3 = DPS48o( oA3, oB3, nmin, peeling, level);
iC4 = DPS48o( oA4, oB4, nmin, peeling, level);
iC5 = DPS48o( oA5, oB5, nmin, peeling, level);
iC6 = DPS48o( oA6, oB6, nmin, peeling, level);
iC7 = DPS48o( oA7, oB7, nmin, peeling, level);
iC8 = DPS48o( oA8, oB8, nmin, peeling, level);
iC9 = DPS48o( oA9, oB9, nmin, peeling, level);
iC10 = DPS48o( oA10, oB10, nmin, peeling, level);
iC11 = DPS48o( oA11, oB11, nmin, peeling, level);
iC12 = DPS48o( oA12, oB12, nmin, peeling, level);
iC13 = DPS48o( oA13, oB13, nmin, peeling, level);
iC14 = DPS48o( oA14, oB14, nmin, peeling, level);
iC15 = DPS48o( oA15, oB15, nmin, peeling, level);
iC16 = DPS48o( oA16, oB16, nmin, peeling, level);
iC17 = DPS48o( oA17, oB17, nmin, peeling, level);
iC18 = DPS48o( oA18, oB18, nmin, peeling, level);
iC19 = DPS48o( oA19, oB19, nmin, peeling, level);
iC20 = DPS48o( oA20, oB20, nmin, peeling, level);
iC21 = DPS48o( oA21, oB21, nmin, peeling, level);
iC22 = DPS48o( oA22, oB22, nmin, peeling, level);
iC23 = DPS48o( oA23, oB23, nmin, peeling, level);
iC24 = DPS48o( oA24, oB24, nmin, peeling, level);
iC25 = DPS48o( oA25, oB25, nmin, peeling, level);
iC26 = DPS48o( oA26, oB26, nmin, peeling, level);
iC27 = DPS48o( oA27, oB27, nmin, peeling, level);
iC28 = DPS48o( oA28, oB28, nmin, peeling, level);
iC29 = DPS48o( oA29, oB29, nmin, peeling, level);
iC30 = DPS48o( oA30, oB30, nmin, peeling, level);
iC31 = DPS48o( oA31, oB31, nmin, peeling, level);
iC32 = DPS48o( oA32, oB32, nmin, peeling, level);
iC33 = DPS48o( oA33, oB33, nmin, peeling, level);
iC34 = DPS48o( oA34, oB34, nmin, peeling, level);
iC35 = DPS48o( oA35, oB35, nmin, peeling, level);
iC36 = DPS48o( oA36, oB36, nmin, peeling, level);
iC37 = DPS48o( oA37, oB37, nmin, peeling, level);
iC38 = DPS48o( oA38, oB38, nmin, peeling, level);
iC39 = DPS48o( oA39, oB39, nmin, peeling, level);
iC40 = DPS48o( oA40, oB40, nmin, peeling, level);
iC41 = DPS48o( oA41, oB41, nmin, peeling, level);
iC42 = DPS48o( oA42, oB42, nmin, peeling, level);
iC43 = DPS48o( oA43, oB43, nmin, peeling, level);
iC44 = DPS48o( oA44, oB44, nmin, peeling, level);
iC45 = DPS48o( oA45, oB45, nmin, peeling, level);
iC46 = DPS48o( oA46, oB46, nmin, peeling, level);
iC47 = DPS48o( oA47, oB47, nmin, peeling, level);
z55 = iC15-iC24;
z49 = iC44+iC15+iC24;
z54 = iC5-iC25;
z48 = iC47-iC25-iC5;
z53 = iC9+iC7;
z51 = iC46+iC9-iC7;
z50 = iC28+iC2;
z52 = iC35+iC28-iC2;
z29 = iC33-iC16;
z47 = iC33+iC19;
z38 = iC32+iC10;
z46 = iC32+iC20;
z28 = iC36+iC18;
z39 = iC36-iC31;
z45 = iC39+iC3;
z26 = iC39+iC22;
z43 = iC40-iC0;
z36 = iC40+iC17;
z34 = iC42+iC14;
z41 = iC42+iC27;
z31 = iC43+iC29;
z40 = iC43+iC30;
z33 = iC45-iC1;
z30 = iC45+iC21;
z19 = z29-z36;
z12 = z29+z36;
z42 = z54+z53;
z25 = z54-z53;
z35 = z55+z50;
z44 = z55-z50;
z32 = z52+z49;
z37 = z52-z49;
z27 = z51+z48;
z63 = z51-z48;
z24 = z39+z38;
z59 = z39-z38;
z13 = z47+z43;
z23 = z47-z43;
z14 = z31+z26;
z22 = z31-z26;
z20 = z41+z30;
z60 = z41-z30;
z11 = z45+z40;
z10 = z45-z40;
z62 = z37+z25;
z61 = z37-z25;
z17 = z34+z33;
z15 = z34-z33;
z57 = z35+z27;
z16 = z35-z27;
z56 = z46+z28;
z58 = z46-z28;
z18 = z42+z32;
z21 = z42-z32;
z64 = z63-z44;
z65 = z63+z44;
oC1 = (z23-z20-z62)/4;
oC3 = (z15-z13-z57)/4;
oC5 = (z60-z57-z12)/4;
oC7 = (z19-z17-z62)/4;
oC8 = (z24-z16+z14)/4;
oC10 = (z61+z58-z22)/4;
oC12 = (z10+z61-z59)/4;
oC14 = (z11+z56-z16)/4;
oC2 = (z18-z58-z23)/4+iC26+iC34;
oC9 = (z20+z22+z21)/4-iC11-iC34;
oC15 = (z17-z10-z21)/4+iC8-iC41;
oC4 = (z59-z19-z18)/4-iC13+iC41;
oC11 = (z15+z14-z64)/4+iC4+iC37;
oC0 = (z24+z65-z13)/4-iC12-iC37;
oC6 = (z56-z65-z12)/4+iC23+iC38;
oC13 = (z11+z60+z64)/4-iC6-iC38;

C = [ oC0 oC1 oC2 oC3 ; oC4 oC5 oC6 oC7 ; oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 ] ;
  end
end
end
