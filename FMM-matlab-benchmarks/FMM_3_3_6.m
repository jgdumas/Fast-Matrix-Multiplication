function C = FMM_3_3_6(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMM_3_3_6(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;3;6>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = A(r2,c0)+A(r0,c0);
t10 = A(r2,c1)-A(r1,c1);
t11 = A(r0,c0)-A(r2,c0);
t12 = A(r0,c1)+A(r2,c1);
t13 = -A(r0,c2)-t9;
t14 = A(r2,c2)-A(r0,c2);
oA10 = A(r1,c2)+t10-t13;
oA16 = t11-t14;
oA17 = t13-A(r2,c2);
oA19 = t14+t11;
oA20 = A(r0,c1)-A(r2,c1)+t11;
oA21 = t9-t12;
oA22 = t9+t12;
oA27 = A(r2,c0)+t10-A(r1,c0);
oA30 = A(r2,c2)+t10-A(r1,c2);
oA2 = -oA17-oA30;
v41 = oA20-oA16;
oA32 = oA2-oA21;
oA31 = -oA10-oA17;
oA34 = oA22-oA10;
oA14 = oA30-oA19;
oA23 = v41-oA19;
oA35 = -oA20-oA14;
oA37 = oA27-oA10;
oA33 = oA31-v41;
oA7 = oA32-oA22;
v51 = -oA31-oA16;
oA12 = oA32-oA23;
oA11 = oA32-oA20;
oA39 = v51-oA27;
oA24 = oA20-oA22+oA27;
oA5 = oA35-oA21;
oA8 = oA33-oA21;
oA1 = -v51;
oA28 = v41-oA34;
oA26 = oA37+oA2;
oA4 = oA34-oA20;
oA36 = oA14-oA27;
oA9 = -oA16-oA30;
oA25 = oA37+oA12;
oA18 = oA17+oA21+oA22;
oA15 = -oA21-oA34;
oA6 = oA31-oA19;
oA3 = oA22+oA35;
oA29 = -oA11-oA19;
oA0 = oA34-oA23;
oA38 = -oA7-oA27;
oA13 = oA22+oA33;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t20 = B(r0,c0)+B(r2,c1)/8;
t21 = B(r1,c0)-B(r2,c3)/8;
r22 = (B(r0,c5)-B(r1,c1))/8;
t22 = B(r2,c4)+r22;
r23 = B(r0,c1)/8;
t24 = B(r2,c1)+B(r2,c3);
t25 = B(r2,c4)-r22;
r26 = B(r1,c3)/8;
t26 = t20-r26;
r27 = (B(r0,c3)+B(r1,c2))/8;
t27 = B(r2,c0)-r27;
r28 = B(r2,c2)/8;
t28 = B(r0,c4)+r28;
r29 = B(r1,c5)/8;
t29 = B(r2,c4)+r23-r29;
t30 = B(r0,c0)+r28;
t31 = B(r2,c2)+B(r2,c5);
t32 = B(r1,c0)+B(r2,c5)/8;
t33 = B(r1,c0)+r29;
t34 = t28-B(r1,c2)/8;
t36 = B(r1,c4)+(B(r0,c5)-B(r2,c5))/8;
t37 = t20+r26;
r39 = B(r0,c2)/8;
t39 = t37+r39;
t40 = t21-B(r0,c3)/8-t36;
t41 = t26-r39;
t43 = t21+t29;
r44 = B(r1,c1)/8;
t44 = t33-r44;
t45 = t22-t27;
t47 = -t25-t27;
t50 = B(r0,c4)-B(r0,c0)-r39+r23;
t51 = r27-t30;
t52 = t21-t29;
t53 = B(r1,c4)-t34-t26;
t55 = t30+r27;
r56 = (t24+t31)/8;
r57 = (t24-t31)/8;
oB4 = t25-t32-t41;
oB5 = t22+t32-t41;
oB6 = t32-t39-t22;
oB7 = t25+t32+t39;
oB8 = t52+t51;
oB9 = t43-t55;
oB10 = t51-t52;
oB11 = t43+t55;
oB14 = (B(r1,c3)-B(r0,c1)-B(r1,c2))/8-B(r2,c0)-t28-t36;
oB16 = r56+t47;
oB18 = t47-r56;
oB21 = t44+t53;
oB26 = t50+t40;
oB27 = t40-t50;
oB29 = t45-r57;
oB31 = t45+r57;
oB33 = r44-B(r1,c4)+t33-t34-t37;
oB34 = t53-t44;
v40 = oB31-oB10;
oB15 = v40+oB34;
v42 = oB7-oB11;
v43 = oB33-oB6;
v44 = oB5-oB9;
oB1 = oB31-v43;
oB2 = oB26+oB15-oB21;
oB12 = v44+oB2;
oB3 = v42-oB14;
v50 = oB8-oB6;
v51 = oB4+v40;
oB20 = oB14-oB1-oB27;
v54 = oB4+oB16;
oB39 = oB9+oB16-oB1;
oB24 = oB4-oB11-oB20;
oB22 = oB10+oB27-oB7;
oB35 = oB29+oB3-oB5;
oB32 = oB11-oB29+oB12;
oB23 = oB26-oB6-oB9;
oB17 = oB10+v44-v54;
oB19 = v50-oB18-v42;
oB13 = v43-v51;
oB25 = oB5+oB8-oB21;
oB0 = oB15-v50;
oB37 = oB12-v54;
oB28 = v50-v51;
oB36 = oB18+oB3-oB8;
oB38 = oB7+oB18-oB15;
oB30 = v44-oB29-v42;

iC0 = FMM( oA0, oB0, nmin, peeling, level);
iC1 = FMM( oA1, oB1, nmin, peeling, level);
iC2 = FMM( oA2, oB2, nmin, peeling, level);
iC3 = FMM( oA3, oB3, nmin, peeling, level);
iC4 = FMM( oA4, oB4, nmin, peeling, level);
iC5 = FMM( oA5, oB5, nmin, peeling, level);
iC6 = FMM( oA6, oB6, nmin, peeling, level);
iC7 = FMM( oA7, oB7, nmin, peeling, level);
iC8 = FMM( oA8, oB8, nmin, peeling, level);
iC9 = FMM( oA9, oB9, nmin, peeling, level);
iC10 = FMM( oA10, oB10, nmin, peeling, level);
iC11 = FMM( oA11, oB11, nmin, peeling, level);
iC12 = FMM( oA12, oB12, nmin, peeling, level);
iC13 = FMM( oA13, oB13, nmin, peeling, level);
iC14 = FMM( oA14, oB14, nmin, peeling, level);
iC15 = FMM( oA15, oB15, nmin, peeling, level);
iC16 = FMM( oA16, oB16, nmin, peeling, level);
iC17 = FMM( oA17, oB17, nmin, peeling, level);
iC18 = FMM( oA18, oB18, nmin, peeling, level);
iC19 = FMM( oA19, oB19, nmin, peeling, level);
iC20 = FMM( oA20, oB20, nmin, peeling, level);
iC21 = FMM( oA21, oB21, nmin, peeling, level);
iC22 = FMM( oA22, oB22, nmin, peeling, level);
iC23 = FMM( oA23, oB23, nmin, peeling, level);
iC24 = FMM( oA24, oB24, nmin, peeling, level);
iC25 = FMM( oA25, oB25, nmin, peeling, level);
iC26 = FMM( oA26, oB26, nmin, peeling, level);
iC27 = FMM( oA27, oB27, nmin, peeling, level);
iC28 = FMM( oA28, oB28, nmin, peeling, level);
iC29 = FMM( oA29, oB29, nmin, peeling, level);
iC30 = FMM( oA30, oB30, nmin, peeling, level);
iC31 = FMM( oA31, oB31, nmin, peeling, level);
iC32 = FMM( oA32, oB32, nmin, peeling, level);
iC33 = FMM( oA33, oB33, nmin, peeling, level);
iC34 = FMM( oA34, oB34, nmin, peeling, level);
iC35 = FMM( oA35, oB35, nmin, peeling, level);
iC36 = FMM( oA36, oB36, nmin, peeling, level);
iC37 = FMM( oA37, oB37, nmin, peeling, level);
iC38 = FMM( oA38, oB38, nmin, peeling, level);
iC39 = FMM( oA39, oB39, nmin, peeling, level);
v49 = iC39+iC35;
v47 = iC0+iC19;
b11 = iC1+iC20-iC27;
v50 = iC34+iC38;
b20 = iC2-iC32;
b21 = iC17+iC30-b20;
b23 = iC11+iC38+v47;
b25 = iC30-iC19-b11;
b29 = iC35-iC22-iC3-b23;
v40 = iC17-b29;
z37 = iC37+iC32+v50;
z36 = iC36+iC33+v49;
z31 = iC31+iC26-iC27-iC30;
z29 = iC29+iC25+b23;
z24 = iC24+iC28-b23;
z18 = iC18-iC21-iC17-iC22;
z16 = iC16+iC23+iC20+iC19;
z15 = iC15+iC21-iC26-iC38+b20;
z14 = iC14+b11-v49;
z13 = iC13+iC39+iC28+b29;
z12 = iC12+iC23-iC0-iC25-iC32-iC38;
z10 = iC10+iC27+v50+b21;
z9 = iC9+iC23-iC39-b25;
z8 = iC8+iC21-iC39-iC33+v40;
z7 = iC7-iC22-b21;
z6 = iC6+b25-iC26-iC33;
z5 = iC5-iC25-v40;
z4 = iC4+iC20-iC34-iC28+v47;
t50 = z37+z24;
b36 = z36-z37+z24;
t45 = z29+z18;
t44 = z29-z18;
b37 = -z31-z16;
t47 = z31-z16;
t42 = z15-z14;
t39 = z15+z14;
t52 = z13-z12;
t41 = z13+z12;
b38 = -z24-z9;
t48 = z10+z8;
b39 = z8-z10;
b40 = -z7-z6;
t46 = z7-z6;
b41 = z4-z5;
t51 = -z5-z4;
t32 = b41-b40;
t28 = b41+b40;
t37 = z9+t48;
t34 = -t51-t46;
t24 = t47+t45;
t27 = b37-t44;
t23 = b37+t44;
oC3 = z9-b39;
b44 = z37+z36+b39;
b45 = t41-t42;
t31 = t42+t41;
b46 = t50-t52;
b47 = -z36-t39;
oC7 = z36+t50+t32;
oC6 = (b36-t31)/8;
oC5 = b45+t32;
oC12 = (t45-t47-oC3-t32)/8;
b51 = b44-b38;
oC10 = (z36-t50+t37-t28)/8;
oC14 = t28+t27;
oC8 = b38+b44+b45;
b53 = b47-b46;
oC13 = t37+t31+t24;
t22 = -t23;
oC2 = t51-t46-t23;
b56 = t48-z9+t22;
oC16 = (t22-b45)/8;
oC17 = t27+b51;
oC4 = (-t34-b51)/8;
oC15 = t34+t24-b53;
oC0 = b53/8;
oC11 = b36+b56;
oC1 = t52+t39+b56;
oC9 = b47+b46-oC2;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
