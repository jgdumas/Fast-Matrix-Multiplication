function C = FMMa_6_3_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
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
    mu=nmin*2^l; ku=nmin*2^l; nu=nmin*2^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMMa_6_3_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMMa(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMMa(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMMa(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<6;3;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t18 = A(r2,c2)+A(r4,c0);
t19 = A(r0,c0)+A(r1,c2);
t20 = A(r2,c1)+A(r3,c0);
t21 = A(r3,c1)-t18;
t22 = A(r0,c1)-A(r1,c1);
t23 = A(r2,c1)+A(r4,c1);
t24 = A(r5,c0)+A(r5,c2);
t25 = A(r0,c2)+A(r1,c0);
t27 = A(r5,c0)-A(r5,c2);
t28 = A(r1,c0)+A(r4,c2)-A(r5,c1);
t29 = A(r0,c0)+A(r2,c2);
t30 = A(r0,c1)-A(r3,c2);
t31 = t19-t21;
t33 = A(r3,c1)+t18;
t34 = A(r4,c2)+A(r2,c0)+A(r3,c1);
t35 = A(r2,c1)-A(r4,c1);
r56 = t19/4;
r57 = t27/4;
r58 = (t22+t23)/4;
r59 = (t28+t29)/4;
r60 = (A(r2,c2)+A(r1,c2)+A(r3,c2))/4;
r61 = (A(r5,c1)-t31)/4;
r62 = (A(r0,c1)+A(r1,c1))/4;
r63 = (t25-t33)/4;
r64 = (t20-t30)/4;
r65 = (t21-t25)/4;
r66 = (t28-t29)/4;
r67 = -(t22+t24)/4;
r68 = (t19-t34)/4;
r69 = -(A(r5,c1)+t31)/4;
r70 = (t24+t35)/4;
r71 = (t23+t27)/4;
r72 = (A(r0,c2)+A(r1,c1)+A(r4,c2)-t20)/4;
r73 = (t22-t23)/4;
r74 = (t20+t30)/4;
oA4 = r67-r68;
oA5 = -r67-r68;
oA7 = t34/4+r56-r57+r62;
oA8 = r74-r59;
oA9 = r66-r64;
oA10 = r64+r66;
oA11 = r59+r74;
oA12 = r70-r63;
oA13 = r63+r70;
oA14 = r65-r71;
oA15 = r65+r71;
oA16 = t24/4+r60-r72;
oA18 = r57-r60-r72;
oA20 = r61-r58;
oA21 = r58+r61;
oA23 = r56-r62+(A(r5,c1)+t33+t35)/4;
oA34 = r69-r73;
oA35 = -r69-r73;
v40 = oA5-oA12;
v41 = v40-oA21;
v42 = oA4+oA13;
v43 = oA7-oA14;
v44 = -v41-oA23;
oA3 = v43-oA11;
oA1 = oA10-v42;
oA0 = v44-oA8;
v48 = oA9+oA16;
oA38 = oA7+oA18-oA15;
oA6 = v44-oA15;
oA29 = oA5+oA35-oA3;
oA22 = v42-v43-oA20;
oA24 = oA4-oA11-oA20;
oA37 = oA12-oA16-oA4;
oA2 = oA9-v40;
oA33 = v42+v44-oA34;
oA27 = oA14-oA20-oA1;
oA30 = -oA9-oA14-oA35;
oA25 = oA5+oA8-oA21;
oA19 = oA11-oA0-oA38;
oA36 = oA18+oA3-oA8;
oA17 = oA5+oA10-v48-oA4;
oA31 = oA10+oA15-oA34;
oA32 = v43-oA35-v40;
oA26 = oA9-oA15-v41;
oA39 = v48-oA1;
oA28 = oA34-oA0-oA4;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t9 = B(r0,c2)+B(r0,c0);
t10 = B(r0,c2)-B(r0,c0);
t11 = B(r1,c2)-B(r1,c1);
t12 = B(r1,c0)-B(r1,c2);
t13 = t11-B(r2,c1);
t14 = B(r2,c0)+t9;
t15 = B(r1,c0)+B(r1,c2);
t16 = t10-B(r2,c0);
oB2 = t14-t13;
oB9 = t16-t13;
oB10 = B(r2,c1)+t11+t14;
oB16 = -t16-B(r2,c2);
oB20 = t12-t10;
oB21 = t9-t15;
oB22 = t9+t15;
oB23 = t10+t12;
oB27 = B(r0,c2)+t11-B(r0,c1);
oB32 = oB2-oB21;
v41 = oB9+oB10;
oB1 = oB2-v41;
oB34 = oB22-oB10;
oB35 = oB9-oB23;
oB33 = oB1-oB20;
oB37 = oB27-oB10;
oB12 = oB32-oB23;
v48 = oB22-oB32;
v49 = oB16+oB9;
v50 = oB20-oB16;
oB14 = -oB35-oB20;
oB11 = oB32-oB20;
oB31 = oB1-oB16;
oB17 = v49-oB2;
oB15 = -oB21-oB34;
oB3 = oB22+oB35;
oB0 = oB34-oB23;
oB38 = v48-oB27;
oB25 = oB37+oB12;
oB19 = v50-oB23;
oB24 = oB20-oB22+oB27;
oB8 = oB11-v41;
oB5 = oB35-oB21;
oB4 = oB34-oB20;
oB7 = -v48;
oB13 = oB22+oB33;
oB18 = v48+v49;
oB28 = v50-oB34;
oB36 = oB14-oB27;
oB39 = -oB27-oB1;
oB6 = oB23+oB33;
oB30 = -v49;
oB26 = oB2+oB37;
oB29 = oB16-oB12;

iC0 = FMMa( oA0, oB0, nmin, peeling, level);
iC1 = FMMa( oA1, oB1, nmin, peeling, level);
iC2 = FMMa( oA2, oB2, nmin, peeling, level);
iC3 = FMMa( oA3, oB3, nmin, peeling, level);
iC4 = FMMa( oA4, oB4, nmin, peeling, level);
iC5 = FMMa( oA5, oB5, nmin, peeling, level);
iC6 = FMMa( oA6, oB6, nmin, peeling, level);
iC7 = FMMa( oA7, oB7, nmin, peeling, level);
iC8 = FMMa( oA8, oB8, nmin, peeling, level);
iC9 = FMMa( oA9, oB9, nmin, peeling, level);
iC10 = FMMa( oA10, oB10, nmin, peeling, level);
iC11 = FMMa( oA11, oB11, nmin, peeling, level);
iC12 = FMMa( oA12, oB12, nmin, peeling, level);
iC13 = FMMa( oA13, oB13, nmin, peeling, level);
iC14 = FMMa( oA14, oB14, nmin, peeling, level);
iC15 = FMMa( oA15, oB15, nmin, peeling, level);
iC16 = FMMa( oA16, oB16, nmin, peeling, level);
iC17 = FMMa( oA17, oB17, nmin, peeling, level);
iC18 = FMMa( oA18, oB18, nmin, peeling, level);
iC19 = FMMa( oA19, oB19, nmin, peeling, level);
iC20 = FMMa( oA20, oB20, nmin, peeling, level);
iC21 = FMMa( oA21, oB21, nmin, peeling, level);
iC22 = FMMa( oA22, oB22, nmin, peeling, level);
iC23 = FMMa( oA23, oB23, nmin, peeling, level);
iC24 = FMMa( oA24, oB24, nmin, peeling, level);
iC25 = FMMa( oA25, oB25, nmin, peeling, level);
iC26 = FMMa( oA26, oB26, nmin, peeling, level);
iC27 = FMMa( oA27, oB27, nmin, peeling, level);
iC28 = FMMa( oA28, oB28, nmin, peeling, level);
iC29 = FMMa( oA29, oB29, nmin, peeling, level);
iC30 = FMMa( oA30, oB30, nmin, peeling, level);
iC31 = FMMa( oA31, oB31, nmin, peeling, level);
iC32 = FMMa( oA32, oB32, nmin, peeling, level);
iC33 = FMMa( oA33, oB33, nmin, peeling, level);
iC34 = FMMa( oA34, oB34, nmin, peeling, level);
iC35 = FMMa( oA35, oB35, nmin, peeling, level);
iC36 = FMMa( oA36, oB36, nmin, peeling, level);
iC37 = FMMa( oA37, oB37, nmin, peeling, level);
iC38 = FMMa( oA38, oB38, nmin, peeling, level);
iC39 = FMMa( oA39, oB39, nmin, peeling, level);

b3 = iC30+iC27;
v49 = iC8+iC18;
b11 = iC32+iC37;
b14 = iC9-iC39-iC30;
b17 = iC38+iC0;
b19 = iC11+iC19+b17;
v47 = iC38+b11;
b21 = iC10-iC22-b11;
b23 = v49-iC33-iC22;
b24 = iC27-iC1-b14;
b25 = iC23-iC33+b14;
b30 = b24-iC23;
z35 = iC35+iC39+iC33+iC36;
z34 = iC34+v47;
z29 = iC29+iC25+b19;
z26 = iC26+iC31-b3;
z24 = iC24+iC28-b19;
z21 = iC21+iC17-iC18+iC22;
z20 = iC20+iC19-b30;
z16 = iC16+b24;
z15 = iC15+iC18+iC31-iC38+b21;
z14 = iC14+iC36-iC19-b25;
z13 = iC13+iC28+b23;
z12 = iC12+iC23-iC32-iC25-b17;
z7 = iC7+iC27+b21;
z6 = iC6+iC31-iC27+b25;
z5 = iC5+b23-iC17-iC39-iC25;
z4 = iC4+iC0-iC28+v47+b30;
z3 = iC3+iC36+v49+b19;
z2 = iC2+iC37-iC17-iC10-b3;
t76 = z24/2;
t75 = z16/2;
t56 = (z3+z2)/2;
t55 = (z29-z16)/2;
t54 = (z13-z12)/2;
t53 = (z29+z16)/2;
t52 = (z3-z2)/2;
t51 = (z6-z7)/2;
t50 = (z14-z15)/2;
t49 = (z5+z4)/2;
t48 = (z35-z21)/2;
t47 = (z34-z20)/2;
t46 = (-z7-z6)/2;
t44 = (z4-z5)/2;
t43 = (-z26-z24)/2;
t42 = (z34+z20)/2;
t40 = (-z13-z12)/2;
t39 = (-z15-z14)/2;
t37 = (-z35-z21)/2;
t36 = (z26-z24)/2;
b42 = t52-z29/2;
b43 = z26/2-t51;
t30 = t52-t43;
oC9 = t48-t56+t42;
t28 = t50-t40;
b45 = t48-t42+t39;
b46 = t46-t44;
b47 = t37-t47;
b48 = t47-t49+t37;
b49 = t56-t36;
t29 = t54-t55;
t24 = t39+b43;
t21 = t44+b47;
b53 = t50+t40+b47;
oC17 = b49-t53;
oC16 = t30-t55;
b54 = t76+t49-t29;
oC7 = t36+t28;
oC14 = t56-t55-t28;
oC4 = b49-b46;
oC2 = t53+b46;
oC12 = b43-t76+b48;
oC6 = t75+t51+b42+b48;
oC3 = t29-b45;
oC0 = t54-b49+b45;
t19 = t46+t21;
oC15 = t21-t46+t28;
oC11 = b54-t24;
oC10 = t24+b54;
oC13 = t43-t19;
oC8 = b42-t75+t19;
oC5 = t55-b53;
oC1 = t30+b53;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
