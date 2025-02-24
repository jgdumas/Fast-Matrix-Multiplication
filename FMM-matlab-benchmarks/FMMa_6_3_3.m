function C = FMMa_6_3_3(A, B, nmin, peeling, level)
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
t18 = A(r2,c1)+A(r3,c0);
t19 = A(r1,c1)-A(r5,c0);
t20 = A(r0,c0)-A(r0,c1);
t21 = A(r1,c2)+A(r5,c2);
t22 = A(r4,c2)-t18;
t23 = A(r2,c2)+A(r3,c2);
t24 = A(r0,c0)+A(r0,c1);
t26 = A(r3,c2)+A(r4,c1);
t27 = A(r3,c1)+A(r2,c0)+A(r4,c2);
t28 = A(r1,c0)-A(r5,c1);
t29 = t19+t22;
t30 = A(r5,c2)-t20;
t31 = A(r1,c2)+A(r4,c0);
t32 = -A(r2,c2)-t24;
r60 = (t19+t27)/4;
r61 = -(t21+t29)/4;
r62 = (A(r0,c2)-t23)/4;
r63 = (A(r4,c1)+A(r1,c1)-A(r2,c1))/4;
r64 = (A(r2,c0)+A(r4,c0)+A(r5,c0))/4;
r65 = (A(r3,c2)+t32)/4;
r66 = (t21-t29)/4;
r67 = -(t22+t28)/4;
r68 = (A(r3,c0)+t26-A(r1,c1))/4;
r69 = (A(r1,c0)-A(r3,c0)+t26-t30)/4;
r70 = (-t19+t27)/4;
r71 = (A(r4,c2)+t18+t28)/4;
r72 = (A(r0,c2)+A(r2,c0)+A(r5,c1)-t31)/4;
r73 = (t21+t24)/4;
r74 = (t20+t23)/4;
r75 = (A(r3,c1)+A(r5,c1)+t31)/4;
r76 = (A(r0,c2)+t23)/4;
r77 = (A(r1,c2)-t30)/4;
oA1 = r68-r72;
oA2 = -r68-r72;
oA4 = r60-r73;
oA5 = r70-r77;
oA6 = -r70-r77;
oA7 = r60+r73;
oA8 = r67-r74;
oA9 = -r67-r74;
oA10 = r65+r71;
oA11 = r71-r65;
oA16 = r66-r62;
oA18 = r61-r76;
oA22 = t32/4-r63-r75;
oA23 = r75+(A(r2,c2)+t20)/4-r63;
oA25 = r64-r69;
oA28 = r62+r66;
oA30 = r76+r61;
oA36 = r64+r69;
v40 = oA5-oA9;
v41 = oA8-oA18;
oA12 = oA2+v40;
v43 = oA10-oA4;
oA3 = v41+oA36;
v45 = oA7-oA11;
oA0 = oA12-oA25-oA23;
oA13 = v43-oA1;
oA14 = v45-oA3;
v50 = oA6-v41;
v51 = oA10-oA22;
v52 = oA8-oA28;
oA31 = v52-oA6+v43;
oA27 = oA7-v51;
oA19 = -v45-v50;
oA20 = v51-oA1-oA11-oA3;
oA21 = oA5+oA8-oA25;
oA32 = oA2+oA7+oA30;
oA35 = -oA14-oA9-oA30;
oA17 = v40+v43-oA16;
oA29 = v40-v45-oA30;
oA39 = oA9+oA16-oA1;
oA38 = oA7-oA0+v50;
oA33 = v52+oA13;
oA34 = oA4+oA28+oA0;
oA24 = oA22+oA3-oA13;
oA37 = oA12-oA4-oA16;
oA15 = oA8+oA0-oA6;
oA26 = oA6+oA9+oA23;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = -B(r1,c2)-B(r1,c1);
t10 = B(r0,c2)-B(r0,c0);
t11 = -B(r0,c2)-B(r0,c0);
t12 = B(r2,c0)+t9;
t13 = B(r2,c1)-t12;
t14 = t10-B(r2,c1);
t15 = B(r1,c0)-B(r1,c2);
oB3 = t13-t11;
oB4 = t14-t12;
oB5 = t13+t11;
oB12 = B(r2,c0)-t9-t14;
oB17 = t11-B(r2,c0)-B(r2,c2);
oB20 = t15-t10;
oB21 = -t11-B(r1,c0)-B(r1,c2);
oB23 = t10+t15;
oB24 = t9-B(r0,c1)-B(r0,c2);
oB32 = oB12+oB23;
oB35 = oB5+oB21;
v42 = oB4+oB12;
oB29 = oB17-oB5;
oB34 = oB20+oB4;
oB22 = oB3-oB35;
oB2 = oB32+oB21;
oB8 = v42-oB3;
v48 = oB20-oB32;
oB10 = oB22-oB34;
oB33 = oB21+oB8;
oB13 = v42-oB5;
oB37 = oB24+oB4;
oB18 = oB3+oB29;
oB28 = -v42-oB29;
oB31 = -oB10-oB17;
oB7 = oB32-oB22;
oB27 = oB24+oB22-oB20;
oB16 = oB12+oB29;
oB1 = oB20+oB33;
oB38 = v48-oB24;
oB9 = oB23+oB35;
oB14 = -oB20-oB35;
oB25 = oB24+v42;
oB11 = -v48;
oB30 = -oB17-oB2;
oB6 = oB4+oB2-oB3;
oB36 = -oB3-oB24;
oB19 = v48-oB29;
oB39 = -oB24-oB13;
oB0 = oB34-oB23;
oB15 = -oB21-oB34;
oB26 = oB2+oB37;

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
b3 = iC39+iC33;
b8 = iC39-iC29-iC5;
v49 = iC28+iC24;
b11 = iC38+iC34;
b12 = iC0+iC4-iC28-iC34;
b13 = iC7-iC30-iC22;
v52 = iC22+iC17;
b18 = iC4+iC24-b11;
v43 = iC33+b8;
b20 = iC26+iC38+b13;
v44 = -b12;
b24 = v44-iC6;
b27 = -iC26-b24;
z37 = iC37+iC32+b11;
z35 = iC35+iC36+b3;
z25 = iC25+iC29+v49;
z23 = iC23+iC16+v44;
z19 = iC19+iC20+b12;
z18 = iC18-iC21-v52;
z15 = iC15+iC21+iC17-b20;
z14 = iC14+iC30+iC36+iC20+b27;
z13 = iC13+iC17-iC24+b8;
z12 = iC12+iC29-iC32-iC16+b18;
z11 = iC11-iC20-b18;
z10 = iC10+iC31+iC34+b20;
z9 = iC9-b3-iC16+b27;
z8 = iC8+iC21+v49-v43;
z3 = iC3+iC36+v52+v43;
z1 = iC1+iC33-iC31+b24;
t76 = (iC27+iC30-iC31-iC26)/2;
t75 = z19/2;
t74 = (iC2+b13-iC32-iC17)/2;
t57 = (z10+z8)/2;
t56 = (-z14-z12)/2;
t55 = (z13-z15)/2;
t54 = (z19+z18)/2;
t53 = (z18-z19)/2;
t52 = (z35+z23)/2;
t51 = (z11-z9)/2;
t50 = (z3-z1)/2;
t49 = (z25-z37)/2;
t47 = (z23-z35)/2;
t43 = (z8-z10)/2;
t42 = (-z11-z9)/2;
t41 = (-z15-z13)/2;
t40 = (z14-z12)/2;
t39 = (z3+z1)/2;
t38 = (z37+z25)/2;
b41 = t55-t56;
t44 = -t52;
b42 = -t74-t50;
t31 = t76+t49;
t46 = -t47;
b45 = t50-t74+t53;
t28 = t51+t43;
b46 = -t44-t42;
b47 = t41-z18/2;
b48 = t49-t76+t40;
b50 = -t74-t39;
b52 = t76+t57-t38;
b53 = t39-t74-t52;
t22 = t57-b47;
t32 = -b42;
oC4 = t76+t38+t32;
t26 = t43-t31;
b54 = t56+t55+t31;
b55 = t40-t75-b46;
oC15 = b41-t46;
oC9 = t52+b42-t28;
oC2 = t54+t28;
oC14 = t54+t32-b41;
t19 = t51-t26;
oC7 = t51+b41+t26;
oC10 = t75+b47-b48;
oC8 = t54+b53;
oC1 = t41+b48+b53;
oC13 = b46+b52;
oC16 = t42+b52+b45;
oC6 = t47+b45;
oC11 = t53-b54;
oC0 = t46+b50+b54;
oC17 = t53-b50-t19;
oC12 = t44+t19;
oC5 = t22-b55;
oC3 = t22+b55;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
