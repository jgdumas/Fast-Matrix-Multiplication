function C = FMM_3_3_6(A, B, nmin, peeling, level)
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
t9 = A(r0,c0)-A(r2,c0);
t10 = A(r0,c2)+A(r1,c2);
t11 = A(r2,c1)+A(r1,c1);
t12 = -A(r2,c0)-A(r0,c0);
t13 = t11-t9;
t15 = t10+A(r2,c1)-A(r1,c1);
t16 = A(r0,c2)-A(r2,c2);
oA4 = t13-t10;
oA10 = t15-t12;
oA11 = t10+t13;
oA15 = t15+t12;
oA16 = t16+t9;
oA17 = t12-A(r0,c2)-A(r2,c2);
oA19 = t9-t16;
oA20 = A(r0,c1)-A(r2,c1)+t9;
oA25 = t11-A(r1,c0)-A(r2,c0);
oA29 = -oA11-oA19;
v41 = oA4+oA15;
v42 = oA17+oA10;
v43 = oA16-oA25;
oA5 = oA17-oA29;
oA37 = oA29-v43;
oA1 = oA16-v42;
oA2 = oA11-v41;
oA28 = -oA4-oA16;
v49 = oA4+oA10;
oA14 = v41-oA5;
v51 = oA19+v43;
oA18 = v42-oA15;
v53 = oA2+oA17;
oA6 = -v42-oA19;
oA12 = oA16-oA29;
oA9 = v53-oA16;
oA3 = oA18-oA29;
oA7 = oA11-v49;
oA36 = v43+v41-v42;
oA38 = oA4+v51;
oA8 = oA1+v41;
oA23 = oA20-oA16-oA19;
oA31 = -v42;
oA30 = -v53;
oA27 = oA10+oA37;
oA35 = -oA20-oA14;
oA26 = -v41-v51;
oA24 = oA37-oA4;
oA32 = oA11+oA20;
oA34 = oA4+oA20;
oA13 = -oA28-oA17;
oA39 = oA5-oA25;
oA22 = oA20+v49;
oA21 = -v41-oA20;
oA33 = oA1-oA20;
oA0 = oA19-oA28;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = B(r0,c3)+B(r1,c2);
r19 = B(r2,c3)/8;
t19 = B(r1,c0)-r19;
t20 = B(r0,c5)-B(r1,c1);
t21 = B(r0,c0)+B(r2,c2)/8;
t22 = t18+t20;
r25 = B(r0,c5)/8;
t25 = B(r0,c4)+r25;
t26 = B(r2,c0)-t22/8;
r27 = B(r1,c3)/8;
t27 = B(r1,c4)-r27;
r28 = B(r0,c1)/8;
r29 = B(r2,c1)/8;
t29 = B(r2,c4)-r29;
r30 = B(r1,c5)/8;
t30 = B(r2,c4)+r28-r30;
r31 = (B(r0,c2)+B(r2,c5))/8;
t31 = B(r0,c0)+r31;
t32 = t19-B(r1,c4);
t33 = B(r2,c4)+r29;
r34 = t18/8;
t34 = t30+r34;
r36 = B(r1,c2)/8;
t36 = t27-r36;
t38 = -t27-r36;
t39 = t31+t32;
t40 = t19+t21;
t42 = B(r2,c0)+r28;
t44 = B(r0,c4)+t21-B(r1,c1)/8+r29;
t45 = r19-t29;
r46 = (B(r2,c2)+B(r2,c5))/8;
t46 = -t26-r46;
t47 = B(r1,c0)+r27;
t48 = t26-r46;
t49 = B(r0,c0)-r31;
t50 = B(r1,c0)-r30;
t51 = t30-r34;
t52 = t33+r19;
t53 = t19-t21;
r54 = t20/8;
r55 = (B(r2,c2)-B(r2,c5))/8;
r56 = (B(r0,c1)+B(r0,c3))/8;
oB5 = t29+t47-t49+r54;
oB7 = t31+t33+t47-r54;
oB8 = t53-t51;
oB9 = t51+t53;
oB10 = t34-t40;
oB11 = t34+t40;
oB12 = r46+t25-t36-t42;
oB14 = t38-t42-r55-t25;
oB16 = t45-t48;
oB17 = B(r2,c0)+t29-r55+(B(r2,c3)+t22)/8;
oB18 = t46-t52;
oB20 = t38-t44-t50;
oB23 = t44-t50-t36;
oB25 = t25+t32-t49+(B(r0,c3)-B(r0,c1))/8;
oB27 = t39-r56-t25;
oB29 = -t45-t48;
oB31 = t52+t46;
oB38 = B(r0,c4)-r25+t39+r56;
v40 = oB7+oB18;
v41 = oB5-oB9;
oB15 = v40-oB38;
v43 = oB12-oB25;
v44 = v41-oB17;
oB0 = v43-oB23;
v47 = oB16-v44;
oB1 = oB14-oB27-oB20;
v49 = oB8-oB15;
v50 = oB11-oB29;
v51 = oB31-oB15;
v52 = oB7-v50;
oB4 = oB10-v47;
oB6 = oB0+v49;
v55 = oB11+oB14;
oB32 = oB12+v50;
oB21 = oB5+oB8-oB25;
oB33 = oB31-oB1+oB6;
oB19 = oB11-oB38-oB0;
oB22 = oB10+oB27-oB7;
oB36 = v40-v55-oB8;
oB13 = v47-oB1;
oB26 = oB9+v43+v49;
oB37 = oB12-v44-oB10;
oB35 = v52-oB5-oB14;
oB24 = oB4-oB11-oB20;
oB30 = v41-v52;
oB39 = oB9+oB16-oB1;
oB2 = oB12-v41;
oB3 = oB7-v55;
oB28 = v47-v51-oB0;
oB34 = oB10-v51;

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

b1 = iC12-iC16;
b2 = iC16+iC20;
b4 = iC31+iC10;
b5 = iC18+iC15;
b6 = iC17-iC5;
v57 = iC18-iC22;
b10 = iC14+iC30;
v51 = iC26-iC30;
b11 = iC25+iC28;
b21 = b4+b5;
b24 = iC13+b6+b11;
b25 = b21-iC37;
b28 = iC20-b1+b25;
b29 = b24-iC35;
b30 = b29-iC26+b10;
b31 = iC11+iC38-b28;
b35 = iC22+b31;
z39 = iC39+b24;
z36 = iC36+iC33-b29;
z34 = iC34+b21-iC22;
z32 = iC32+iC38+iC22-b25;
z29 = iC29+b35;
z27 = iC27-iC31-v51;
z24 = iC24+b11-b35;
z23 = iC23+iC19+b2;
z21 = iC21+iC17-v57;
z9 = iC9+iC35-b2-b10;
z8 = iC8+iC13-iC33+iC28+v57;
z7 = iC7+iC38+v51-b5;
z6 = iC6+b30-iC33-iC19;
z4 = iC4+iC37-b11+b1;
z3 = iC3+iC22*2-iC35+b6+b31;
z2 = iC2+iC37-iC26-iC17-b4;
z1 = iC1+iC20-iC31+b30;
z0 = iC0+iC25-iC22+iC19+b28;
t53 = z39+z36;
b37 = z36-z39;
t58 = z29/8;
b38 = z24-z27;
t40 = z27+z24;
t50 = z34-z23;
t41 = z34+z23;
t49 = z32-z21;
t35 = z32+z21;
t57 = z9+z8;
t55 = z9-z8;
b39 = z29+z7;
b40 = z4-z6;
b41 = z7+z6+z4;
t52 = z3+z1;
b43 = z1-z3;
b44 = z0-z2;
b45 = z2+z0;
b47 = z4-z7+b38;
t38 = -t55;
b48 = b40+b39;
t43 = -b43;
t27 = t49-t41;
b49 = t49-b45+t41;
t29 = t53-t40;
b50 = t53+t43+t40;
t32 = z6+b37;
t28 = b38+b37;
oC16 = (t43+b45)/8-t58;
b51 = t38-z29;
t30 = t50+t35;
oC12 = t58+(t38-b41)/8;
oC5 = b41-t30;
oC1 = b51+t30;
oC8 = t38+t29;
b54 = t52-b44+t28;
t26 = -t27;
b55 = t52+b44+t27;
t23 = t57+t26;
t20 = t32+b47;
oC9 = b48-t29;
oC13 = z29+t23;
oC10 = (t32-b47+t23)/8;
oC11 = b51+b50-b45;
oC7 = b41+b45+b50;
oC3 = t55+b43+b49;
oC2 = t43-b48+b49;
oC15 = z29+t20;
oC4 = (t50-t57-t35-t20)/8;
oC17 = t57-z29+b54;
oC0 = (-t26-b54)/8;
oC14 = b40+b55-b39;
oC6 = (t28+b55)/8;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
