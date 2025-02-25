function C = FMMa_3_3_6(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMMa_3_3_6(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
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
    if l>=level, fprintf("# Core<3;3;6>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = A(r2,c0)-A(r0,c0);
t10 = A(r2,c1)+A(r1,c1);
t11 = A(r0,c2)+A(r1,c2);
t12 = -A(r2,c0)-A(r0,c0);
t13 = t9+t10;
t14 = A(r0,c2)+A(r2,c2);
oA0 = t10-t11-t9;
oA4 = t13-t11;
oA10 = A(r2,c1)+t11-t12-A(r1,c1);
oA11 = t11+t13;
oA17 = t12-t14;
oA18 = -t12-t14;
oA19 = A(r2,c2)-t9-A(r0,c2);
oA23 = A(r0,c1)-A(r2,c1)+t9;
oA25 = t10-A(r1,c0)-A(r2,c0);
v40 = oA11-oA4;
v41 = oA10+oA17;
v42 = oA0-oA19;
oA38 = oA0-oA25;
oA34 = oA23+oA0;
oA7 = v40-oA10;
oA29 = -oA11-oA19;
oA8 = v42-oA18;
oA15 = v41-oA18;
oA6 = -v41-oA19;
oA16 = v42-oA4;
oA5 = oA17-oA29;
oA9 = oA11-oA8-oA10;
oA30 = -oA7-oA18;
oA13 = v42-oA17;
oA26 = -oA38-oA15;
oA3 = oA18-oA29;
oA33 = oA6-oA23;
oA14 = oA30-oA19;
oA22 = oA10+oA34;
oA1 = oA16-v41;
oA27 = -oA38-oA7;
oA32 = oA34+v40;
oA28 = -v42;
oA20 = oA34-oA4;
oA37 = -v40-oA38;
oA36 = oA8-oA25;
oA2 = v40-oA15;
oA35 = oA9-oA23;
oA39 = oA5-oA25;
oA24 = -oA11-oA38;
oA12 = oA0+v40;
oA21 = -oA34-oA15;
oA31 = -v41;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = B(r0,c0)+B(r2,c1);
t19 = B(r0,c4)+B(r2,c2);
t20 = B(r0,c5)+B(r2,c5);
t21 = B(r1,c1)-t20;
t22 = B(r0,c5)-B(r2,c5);
t25 = B(r1,c2)-B(r1,c4);
t26 = B(r1,c0)+B(r1,c1);
t29 = t19+B(r1,c5)+t18;
t30 = B(r0,c2)+B(r1,c3)+B(r2,c4);
t31 = B(r2,c1)-B(r2,c4);
t32 = B(r1,c2)+B(r1,c4);
t33 = B(r2,c0)-B(r2,c2);
t34 = B(r0,c3)+B(r1,c2)-t21;
t35 = B(r1,c3)-B(r0,c1)-B(r2,c0);
r56 = (-t31+t34)/4;
r57 = (-t18+t30)/4;
r58 = (B(r2,c3)-t33)/4;
r59 = (t19+t35)/4;
r60 = (t32+B(r1,c1)-B(r1,c0))/4;
r61 = (t31+t34)/4;
r62 = (B(r0,c3)+B(r0,c4))/4;
r63 = (t20+t25)/4;
r64 = (B(r1,c3)-t29)/4;
r65 = (-t19+t35)/4;
r66 = -(t18+t30)/4;
r67 = (t22-t26)/4;
r68 = (t25-t26)/4;
r69 = (B(r1,c3)+t29)/4;
r70 = (B(r2,c3)+t33)/4;
r71 = (B(r0,c2)-B(r1,c4)-B(r2,c3))/4;
r72 = (B(r1,c0)-t21)/4;
r73 = (t22+t32)/4;
oB1 = r62-r71+(B(r2,c1)-B(r1,c1)-B(r1,c5)-B(r2,c0))/4;
oB4 = r57-r72;
oB5 = r72+r57;
oB6 = r66-r67;
oB7 = -r66-r67;
oB12 = r59+r63;
oB13 = r63-r59;
oB14 = r65-r73;
oB15 = r73+r65;
oB16 = r58+r61;
oB17 = r70+r56;
oB22 = r68-r69;
oB23 = r68+r69;
oB28 = r61-r58;
oB29 = r56-r70;
oB34 = r60+r64;
oB35 = r60-r64;
oB38 = r62+r71+(B(r0,c0)+B(r0,c1)+B(r1,c0)-t22)/4;
v40 = oB34-oB4;
v41 = oB13+oB1;
v42 = oB5+oB35;
oB0 = v40-oB28;
oB3 = v42-oB29;
v46 = v41-oB16-oB17;
v47 = oB7-oB14;
v48 = oB6+oB15;
v50 = oB13-oB22;
v51 = oB5+oB23;
oB10 = oB4+v41;
v53 = oB3+oB38+oB0;
oB26 = oB6+v46+v51;
oB9 = oB5+v46;
oB25 = oB12-oB23-oB0;
oB11 = v47-oB3;
oB2 = oB12+v46;
oB32 = oB12-v42+v47;
oB30 = -v42-v46-oB14;
oB36 = v53-oB6-oB7;
oB27 = oB7+oB22-oB10;
oB21 = v48+v51-oB12;
oB31 = oB15-v40+v41;
oB24 = oB3-v50;
oB33 = oB13-v40+v48;
oB8 = v48-oB0;
oB39 = oB5+oB13-oB17;
oB19 = v47-v53;
oB20 = oB4-v47+v50;
oB37 = oB12-oB16-oB4;
oB18 = oB15+oB38-oB7;

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
b1 = iC12-iC16;
b4 = iC18-iC22;
v47 = iC17-iC5;
v52 = iC30-iC26;
b7 = iC31+iC10;
v54 = iC11+iC29;
b11 = iC16+iC20;
b14 = iC25+iC28;
b17 = iC20-v54-b1;
b25 = iC13+v47+b14;
b26 = b17-iC37;
b27 = b25-iC35;
b28 = iC14+v52+b27;
b29 = b26-iC22+b7;
z39 = iC39+b25;
z36 = iC36+iC33-b27;
z34 = iC34+iC38-b26;
z32 = iC32+b17;
z27 = iC27+v52-iC31;
z24 = iC24+iC29+b14;
z23 = iC23+iC19+b11;
z21 = iC21+iC17-b4;
z9 = iC9+iC35-iC30-iC14-b11;
z8 = iC8+iC13-iC33+iC28+b4;
z7 = iC7+b29-v52;
z6 = iC6+b28-iC33-iC19;
z4 = iC4+iC37-b14+b1;
z3 = iC3+iC22-iC29-iC35+v47;
z2 = iC2+iC37-iC17-iC26-b7;
z1 = iC1+iC20-iC31+b28;
z0 = iC0+iC38+iC19+iC25+v54;
t75 = (iC15+iC18-iC38+b29)/2;
t74 = z7/2;
t56 = (z36-z39)/2;
t55 = (z27-z24)/2;
t54 = (z23-z21)/2;
t53 = (z4-z6)/2;
t52 = (z9+z8)/2;
t51 = (z8-z9)/2;
t50 = (-z3-z1)/2;
t49 = (z34-z32)/2;
t48 = (z34+z32)/2;
t47 = (-z23-z21)/2;
t45 = (-z2-z0)/2;
t44 = (z0-z2)/2;
t42 = (z1-z3)/2;
t41 = (z39+z36)/2;
t39 = (z27+z24)/2;
t33 = (z7+z6+z4)/2;
b43 = z6/2+t56;
t28 = t56-t55;
b44 = (z7-z4)/2+t55;
b45 = t74+t53;
b47 = t48+t75-t54;
b48 = t47-t52;
b49 = t54+t48+t45;
t26 = t41-t39;
t25 = t47-t49;
b52 = t41+t39-t42;
oC16 = t75-t45-t42;
oC12 = t51-t33;
t21 = t49-b48;
b54 = t50+t44-t28;
b56 = b44-b43;
oC9 = b45-t26-t75;
oC8 = t26-t75+t51;
b57 = t44-t50+t25;
oC5 = t33-b47;
oC1 = t51+b47;
oC11 = t51+t45+b52;
oC7 = t33+b52-t45;
oC3 = t42+b49-t51;
oC2 = b49-t42-b45;
oC13 = t21+t75;
oC10 = b43+b44+t21;
oC15 = t75-b56;
oC4 = t49+b48+b56;
oC14 = t53-t74+b57;
oC6 = t28-t75+b57;
oC17 = t52-b54;
oC0 = t25-t75+b54;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
