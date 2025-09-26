function C = DPS48a_mul_4_4_4(A, B, nmin, peeling, level)
NTH2ROOT2once=nthroot(2,2);
NTH2ROOT2o2=nthroot(2,2)/2;
NTH2ROOT2t2o=2/nthroot(2,2);
NTH2ROOT2t1o=1/nthroot(2,2);
NTH2ROOT2o4=nthroot(2,2)/4;
NTH2ROOT2o8=nthroot(2,2)/8;
NTH2ROOT2o16=nthroot(2,2)/16;
NTH2ROOT2o3=nthroot(2,2)/3;
NTH2ROOT2f2o3=nthroot(2,2)*2/3;

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
    C(1:mu,1:nu)=DPS48a_mul_4_4_4(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS48a_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS48a_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS48a_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<4;4;4>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
t16 = A(r2,c0)+A(r2,c2);
t17 = A(r1,c1)+A(r1,c2);
t18 = A(r1,c0)-A(r3,c0);
t19 = A(r0,c0)-A(r1,c0);
r20 = A(r2,c1)*2;
t20 = t17-r20;
t21 = A(r0,c0)-t18;
t22 = A(r0,c3)-A(r3,c0);
t23 = A(r1,c2)-A(r2,c1);
t24 = A(r0,c3)+t21;
r25 = A(r1,c1)*2;
t25 = t16-r25;
r26 = t23*2;
t27 = A(r0,c0)+t18;
t28 = t18+A(r1,c1)+r26;
t29 = A(r0,c3)+t16;
t30 = t16+t22;
t31 = t20+A(r0,c2);
t32 = A(r0,c3)-A(r0,c1);
t33 = t17-A(r2,c1);
t34 = A(r3,c3)-t17;
t35 = A(r1,c2)+t18;
t36 = A(r0,c0)+A(r1,c0);
t37 = t16-t22;
t38 = t24-t25;
r39 = A(r1,c0)*2;
oA0 = A(r2,c0)+r20;
oA2 = A(r2,c3)-A(r3,c0)+t20-t29;
oA4 = A(r0,c3)+A(r3,c1)-t21;
oA6 = t38-A(r3,c3);
oA7 = A(r3,c0)*2-A(r2,c3);
oA8 = A(r1,c0)-A(r1,c1)+A(r2,c1)+t16;
oA9 = t31-t19;
oA10 = t19-t34;
oA11 = -t28;
oA12 = t31-t37;
oA13 = A(r0,c3)*2-A(r0,c1);
oA14 = (t29-t27)/2-t20;
oA15 = A(r0,c2)+t28*2;
oA17 = r25-A(r3,c2);
oA19 = t36-t34;
oA22 = A(r3,c2)-t17-t30-r39;
oA23 = t16+r39-A(r1,c3);
oA24 = t23-A(r0,c3);
oA28 = A(r1,c3)+t25+r20;
oA30 = A(r1,c2)*2+A(r1,c3);
oA32 = A(r3,c1)+r26-t24;
oA33 = t16+t35;
oA34 = A(r0,c0)-A(r3,c0)+t33;
oA36 = A(r1,c1)-A(r1,c2)-A(r3,c2)+t19;
oA37 = (A(r0,c0)-t30)/2-A(r1,c0)*3/2-A(r1,c2);
oA41 = A(r2,c2)+(A(r2,c0)+t35)*2;
oA42 = A(r2,c3)-t20-t36;
oA43 = (t19-t37)/2;
oA44 = t27+t32+t33*2;
oA45 = t19-t32-A(r3,c0);
oA46 = A(r3,c1)-A(r0,c0)*2;
oA47 = t38/2;
oA1 = A(r1,c3);
oA3 = A(r0,c3);
oA5 = A(r3,c2);
oA16 = A(r0,c0);
oA18 = A(r1,c1);
oA20 = A(r2,c2);
oA21 = A(r3,c0);
oA25 = A(r2,c0);
oA26 = A(r3,c1);
oA27 = A(r2,c1);
oA29 = A(r0,c1);
oA31 = A(r1,c0);
oA35 = A(r2,c3);
oA38 = A(r3,c3);
oA39 = A(r1,c2);
oA40 = A(r0,c2);

[m,n] = size(B);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
t16 = B(r0,c3)-B(r0,c2);
oB16 = B(r2,c3)-t16;
oB21 = t16-B(r3,c1);
v48 = B(r2,c0)+B(r2,c1);
v49 = B(r0,c3)-oB21;
v50 = oB16-B(r2,c3);
v51 = B(r3,c2)+B(r2,c2);
oB24 = v48-B(r0,c1);
oB39 = B(r2,c3)-v49;
oB14 = v51-B(r0,c1);
v55 = v50+oB21;
v56 = B(r0,c1)-B(r1,c0);
v57 = B(r3,c2)+B(r2,c2)*2;
oB31 = oB16-v49;
oB11 = B(r0,c3)+v50;
v60 = B(r1,c3)-oB31;
v61 = B(r2,c2)-B(r2,c0);
v62 = B(r3,c0)+oB21;
v63 = oB16-B(r0,c3);
oB37 = v48-B(r2,c2);
oB43 = v51-oB24;
v66 = oB39+B(r1,c1);
v67 = oB14+B(r1,c2);
v68 = B(r3,c3)-oB24;
v69 = v55-B(r0,c0);
oB34 = v57-v48;
oB4 = v60+oB24;
oB26 = B(r1,c1)-v49-v50;
oB41 = v61-v55;
oB36 = oB37+v66;
oB33 = -v55;
oB9 = B(r1,c2)-v48+v51*2;
oB5 = B(r2,c3)-v68;
oB25 = B(r3,c2)-v62;
oB7 = B(r1,c3)-oB21*2-v63;
oB1 = oB39+v67;
oB40 = v56+oB11;
oB17 = B(r2,c3)-v56;
oB38 = B(r2,c2)+v66;
oB45 = v60+oB34;
oB29 = B(r0,c3)+v61;
oB15 = oB11-v68;
oB12 = v69-oB43;
oB44 = B(r0,c1)+B(r3,c3)-v57;
oB19 = v63-B(r3,c0);
oB46 = oB16-v67;
oB2 = oB14-v69;
oB13 = B(r0,c0)+v49-v50;
oB0 = v56-B(r3,c2);
oB23 = B(r0,c1)+v62;
oB3 = B(r0,c3);
oB6 = B(r2,c0);
oB8 = B(r0,c1);
oB10 = B(r1,c3);
oB18 = B(r2,c3);
oB20 = B(r0,c0);
oB22 = B(r2,c1);
oB27 = B(r3,c2);
oB28 = B(r1,c0);
oB30 = B(r1,c1);
oB32 = B(r3,c3);
oB35 = B(r3,c0);
oB42 = B(r1,c2);
oB47 = B(r2,c2);

iC0 = DPS48a_mul( oA0, oB0, nmin, peeling, level);
iC1 = DPS48a_mul( oA1, oB1, nmin, peeling, level);
iC2 = DPS48a_mul( oA2, oB2, nmin, peeling, level);
iC3 = DPS48a_mul( oA3, oB3, nmin, peeling, level);
iC4 = DPS48a_mul( oA4, oB4, nmin, peeling, level);
iC5 = DPS48a_mul( oA5, oB5, nmin, peeling, level);
iC6 = DPS48a_mul( oA6, oB6, nmin, peeling, level);
iC7 = DPS48a_mul( oA7, oB7, nmin, peeling, level);
iC8 = DPS48a_mul( oA8, oB8, nmin, peeling, level);
iC9 = DPS48a_mul( oA9, oB9, nmin, peeling, level);
iC10 = DPS48a_mul( oA10, oB10, nmin, peeling, level);
iC11 = DPS48a_mul( oA11, oB11, nmin, peeling, level);
iC12 = DPS48a_mul( oA12, oB12, nmin, peeling, level);
iC13 = DPS48a_mul( oA13, oB13, nmin, peeling, level);
iC14 = DPS48a_mul( oA14, oB14, nmin, peeling, level);
iC15 = DPS48a_mul( oA15, oB15, nmin, peeling, level);
iC16 = DPS48a_mul( oA16, oB16, nmin, peeling, level);
iC17 = DPS48a_mul( oA17, oB17, nmin, peeling, level);
iC18 = DPS48a_mul( oA18, oB18, nmin, peeling, level);
iC19 = DPS48a_mul( oA19, oB19, nmin, peeling, level);
iC20 = DPS48a_mul( oA20, oB20, nmin, peeling, level);
iC21 = DPS48a_mul( oA21, oB21, nmin, peeling, level);
iC22 = DPS48a_mul( oA22, oB22, nmin, peeling, level);
iC23 = DPS48a_mul( oA23, oB23, nmin, peeling, level);
iC24 = DPS48a_mul( oA24, oB24, nmin, peeling, level);
iC25 = DPS48a_mul( oA25, oB25, nmin, peeling, level);
iC26 = DPS48a_mul( oA26, oB26, nmin, peeling, level);
iC27 = DPS48a_mul( oA27, oB27, nmin, peeling, level);
iC28 = DPS48a_mul( oA28, oB28, nmin, peeling, level);
iC29 = DPS48a_mul( oA29, oB29, nmin, peeling, level);
iC30 = DPS48a_mul( oA30, oB30, nmin, peeling, level);
iC31 = DPS48a_mul( oA31, oB31, nmin, peeling, level);
iC32 = DPS48a_mul( oA32, oB32, nmin, peeling, level);
iC33 = DPS48a_mul( oA33, oB33, nmin, peeling, level);
iC34 = DPS48a_mul( oA34, oB34, nmin, peeling, level);
iC35 = DPS48a_mul( oA35, oB35, nmin, peeling, level);
iC36 = DPS48a_mul( oA36, oB36, nmin, peeling, level);
iC37 = DPS48a_mul( oA37, oB37, nmin, peeling, level);
iC38 = DPS48a_mul( oA38, oB38, nmin, peeling, level);
iC39 = DPS48a_mul( oA39, oB39, nmin, peeling, level);
iC40 = DPS48a_mul( oA40, oB40, nmin, peeling, level);
iC41 = DPS48a_mul( oA41, oB41, nmin, peeling, level);
iC42 = DPS48a_mul( oA42, oB42, nmin, peeling, level);
iC43 = DPS48a_mul( oA43, oB43, nmin, peeling, level);
iC44 = DPS48a_mul( oA44, oB44, nmin, peeling, level);
iC45 = DPS48a_mul( oA45, oB45, nmin, peeling, level);
iC46 = DPS48a_mul( oA46, oB46, nmin, peeling, level);
iC47 = DPS48a_mul( oA47, oB47, nmin, peeling, level);
t37 = iC40+iC35;
t32 = iC41-iC30;
b1 = -iC32-iC28;
b3 = iC32/2-iC24;
b5 = iC35/2+iC21;
t33 = iC19-iC17;
t38 = iC22-iC12;
t36 = iC36+iC9;
b12 = iC10-iC20-iC5;
t30 = iC44-iC0;
b18 = -iC43-(iC12+iC9)/2;
t25 = iC37-(iC36+iC22)/2;
t23 = iC31+(iC10-iC19)/2;
oC1 = iC23+iC4+t38;
oC15 = iC45-iC25+t36;
t21 = iC34+(iC44+t36)/2;
oC6 = iC2-iC6+b1;
b22 = b1/2-iC14;
oC12 = iC46-iC13+t33;
oC3 = iC15+iC7+t32;
t20 = b5-b3;
oC13 = iC42+iC38+t30;
t24 = iC47+t30/2;
t26 = -b18;
t22 = -b22;
oC14 = iC1+b12;
t19 = b18+t25;
b32 = b22-iC6+t24;
oC10 = iC16+iC13-iC3+(iC42+iC6-t37-t33)/2+t26+t22;
oC5 = iC39+(iC42+iC30-b12-t26-t25+t24+t22)/2;
b35 = -t23-t21;
oC11 = iC11-iC7-(iC40+iC25+iC4+t32)/2+t21-t20;
oC8 = iC8+iC4+(iC7-iC28+t38)/2+t23+t20;
oC4 = iC18+(iC25-iC17-iC5-iC4)/2+b3+b35;
oC9 = iC27+iC25+(iC7+iC0)/2+b5+b35;
b39 = iC3*2-iC13+t19+b32;
oC7 = iC33+(iC20-iC41+t19-b32)/2;
oC0 = iC29-b39;
oC2 = iC26+t37+b39;

C = [ oC0 oC1 oC2 oC3 ; oC4 oC5 oC6 oC7 ; oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 ] ;
  end
end
end
