function C = FMMa_3_3_6(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
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
    mu=nmin*2^l; ku=nmin*2^l; nu=nmin*2^l;
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
t9 = -A(r1,c1)-A(r1,c2);
t10 = A(r2,c0)+A(r2,c1);
t11 = A(r0,c0)+A(r0,c2);
t12 = A(r2,c2)+A(r2,c1);
t13 = A(r2,c0)+A(r2,c2);
t14 = A(r0,c2)-t9;
oA7 = t14-A(r0,c0)-t10;
oA13 = A(r1,c1)-A(r1,c2)+t10+t11;
oA16 = t11-t13;
oA17 = -t11-t13;
oA20 = A(r0,c0)+A(r0,c1)-t10;
oA27 = t10-A(r1,c0)-A(r1,c1);
oA29 = t9-t12;
oA30 = t12+t9;
oA32 = A(r0,c1)+t14;
oA11 = oA32-oA20;
v41 = oA7+oA13;
oA1 = v41-oA11;
v43 = -oA29-oA30;
oA31 = oA1-oA16;
oA28 = -oA17-oA13;
oA2 = -oA17-oA30;
oA12 = oA16-oA29;
oA19 = -oA11-oA29;
oA10 = -oA17-oA31;
v50 = oA27+oA7;
oA4 = -oA28-oA16;
oA39 = -oA1-oA27;
oA14 = oA11-v43;
oA5 = oA17-oA29;
oA9 = -oA30-oA16;
oA33 = v41-oA32;
oA22 = oA32-oA7;
oA26 = oA27-oA30+oA31;
oA18 = -oA7-oA30;
oA25 = oA5-oA39;
oA37 = oA27-oA10;
oA21 = oA2-oA32;
oA23 = oA32-oA12;
oA24 = v50-oA11;
oA3 = v43-oA7;
oA6 = v41-oA12;
oA35 = v43-oA32;
oA38 = -v50;
oA15 = oA11-oA9-oA13;
oA34 = oA20+oA4;
oA36 = oA14-oA27;
oA0 = oA19-oA28;
oA8 = v41-oA2;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6; 
t18 = B(r0,c4)+B(r2,c1);
t19 = B(r0,c0)+B(r2,c2);
t20 = B(r1,c0)-B(r1,c1);
t21 = B(r1,c4)+B(r2,c3);
t22 = B(r0,c3)-B(r1,c1);
t23 = B(r1,c2)+B(r1,c4);
t24 = B(r0,c1)+B(r0,c4);
t26 = B(r2,c1)-B(r2,c4);
t28 = B(r0,c5)+B(r2,c5);
t29 = B(r0,c5)-B(r2,c5);
t30 = B(r1,c3)-(t18+t19);
t31 = B(r0,c2)+B(r1,c5)+B(r2,c0);
t32 = B(r0,c3)+B(r1,c0);
t33 = B(r1,c2)+B(r2,c3);
r60 = t19/4;
r61 = -(B(r1,c5)-t30)/4;
r62 = (t21-t22)/4;
r63 = (t23+t29)/4;
r64 = (t18+t31)/4;
r65 = (t20-t23)/4;
r66 = -(t20+t28)/4;
r67 = (B(r1,c3)-B(r2,c2)-(B(r2,c0)+t24))/4;
r68 = (t20+t23)/4;
r69 = (t33+t22+t28)/4;
r70 = (t24+B(r0,c0)+B(r0,c2))/4;
r71 = (B(r2,c0)-(B(r2,c2)+t26))/4;
r72 = (t21+t29-t32)/4;
r73 = (B(r0,c1)-(B(r1,c5)-B(r2,c4)))/4;
r74 = (t18-t31)/4;
r75 = (B(r1,c5)+t30)/4;
r76 = (t26-(B(r0,c2)-(B(r0,c0)-B(r1,c3))))/4;
r77 = (t21+t22)/4;
oB0 = r62-r64;
oB1 = r74+r77;
oB2 = r74-r77;
oB3 = r62+r64;
oB4 = r66-r76;
oB5 = -r66-r76;
oB10 = r73+(B(r0,c3)-B(r1,c0)+t33)/4-r60;
oB11 = r60+r73+(B(r1,c2)-B(r2,c3)+t32)/4;
oB14 = r67-r63;
oB15 = r63+r67;
oB16 = r69-r71;
oB17 = r69+r71;
oB20 = r75-r68;
oB21 = r68+r75;
oB34 = r61-r65;
oB35 = -r61-r65;
oB36 = r70+r72;
oB38 = r70-r72;
v40 = oB10-oB17;
v41 = oB11+oB14;
v42 = oB4-v40;
v43 = oB38-v41;
v44 = oB15+v43;
v45 = v42+oB16;
oB8 = v44-oB36;
oB9 = oB5-v45;
v48 = oB10-oB1;
v49 = oB2-oB9;
v50 = oB3-oB35;
v51 = oB11+oB20;
v52 = oB34-oB0;
v53 = oB21-oB8;
oB7 = oB3+v41;
oB18 = v44-oB3;
oB13 = v48-oB4;
oB22 = v48-v51-oB3;
oB26 = oB2-oB15+oB21;
oB24 = oB4-v51;
oB23 = v49+v53-oB0;
oB39 = oB5-v42-oB1;
oB27 = oB14-oB20-oB1;
oB12 = oB2+v45;
oB33 = oB8+v48-v52;
oB31 = oB10+oB15-oB34;
oB30 = -oB14-oB35-oB9;
oB29 = oB5-v50;
oB32 = oB11+v49+v50;
oB6 = oB0-oB36+v43;
oB25 = oB5-v53;
oB37 = oB2-v40;
oB19 = oB11-oB38-oB0;
oB28 = v52-oB4;

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

b2 = iC0-iC28;
v51 = iC24+iC29;
b4 = iC35+iC5+iC29;
b5 = iC26-iC37-iC2;
b6 = iC6-iC26+iC36;
b14 = iC18-iC21;
b15 = iC37+iC34+iC32;
b20 = iC34-iC4-b2;
b25 = iC24+iC36+b4;
v42 = -b20;
b29 = v42-iC16;
z39 = iC39+iC35+iC36+iC33;
z38 = iC38+b15;
z27 = iC27+iC30-iC26-iC31;
z25 = iC25+iC28+v51;
z23 = iC23-b29;
z22 = iC22+iC17-b14;
z20 = iC20+iC19+v42;
z15 = iC15+iC21+iC34-b5;
z14 = iC14+iC30-iC19+b6;
z13 = iC13+iC17-iC33-b25;
z12 = iC12+iC37+iC4-iC16+v51;
z11 = iC11+iC19-iC24+b2-b15;
z10 = iC10+iC17+iC31+b5;
z9 = iC9+iC35+b6+b29;
z8 = iC8+iC28+iC21+b25;
z7 = iC7+iC2-iC30-iC32-b14;
z3 = iC3+b14-b4;
z1 = iC1+iC33-iC6-iC31+b20;
t78 = z23/2;
t77 = z7/2;
t57 = (z3-z1)/2;
t56 = (z38-z25)/2;
t55 = (z10+z8)/2;
t54 = (z13-z14)/2;
t53 = (-z11-z9)/2;
t52 = (z39-z27)/2;
t51 = (z14+z13)/2;
t50 = (z38+z25)/2;
t49 = (z23+z7)/2;
t47 = (z10-z8)/2;
t46 = (-z22-z20)/2;
t45 = (z39+z27)/2;
t44 = (z15+z12)/2;
t41 = (z3+z1)/2;
t39 = (z11-z9)/2;
t38 = (z12-z15)/2;
t36 = (z20-z22)/2;
t34 = (z23-z7)/2;
b39 = t55-z22/2;
b40 = (z20-z23)/2+t53;
b41 = t77+t56-t52;
t40 = -t50;
b45 = t41-t46;
b46 = -t78-t44-t40;
oC11 = t55+t53+t57+t45+t40;
t30 = t47-t39;
t32 = t51-t38;
t29 = t51+t38;
t27 = t57-t36;
t28 = t50+t45;
b47 = t54+t45-b45;
oC14 = b45-t49;
t24 = -t30;
oC5 = t49+t36+t29;
t21 = t77+t28;
b48 = t47+t39+t28;
oC3 = t78+t30-t27;
oC2 = t34+t27;
t20 = b39-b40;
oC1 = b39+t32+b40;
oC16 = t57-t29;
oC12 = t24-t77;
oC8 = t56+t52+t29+t24;
oC9 = b41-t32;
oC15 = t44-t54-t21;
oC7 = t57+t21;
oC13 = t44+t54+t20;
oC10 = b41+t20;
oC17 = t41-b48;
oC4 = t46-t34+b48;
oC6 = b46-b47;
oC0 = b46+b47;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
