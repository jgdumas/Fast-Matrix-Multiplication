function C = FMMa363_mul_3_6_3(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMMa363_mul_3_6_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMMa363_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMMa363_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMMa363_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;6;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = A(r0,c0)-A(r0,c2);
t19 = A(r2,c4)-A(r0,c3);
oA6 = A(r1,c3)-A(r2,c4)-A(r1,c0);
t22 = A(r1,c4)+A(r2,c2);
oA10 = A(r0,c5)+t18;
oA3 = t19-A(r1,c5);
t25 = A(r0,c1)-A(r0,c4);
t26 = A(r2,c0)-t22;
t27 = A(r2,c5)+A(r2,c3);
oA0 = t22-A(r2,c5);
oA9 = oA6-t25;
oA16 = t25+t26;
oA18 = t18+t19-A(r1,c2);
oA20 = A(r0,c1)-A(r1,c4)+A(r2,c0);
oA21 = A(r0,c2)+A(r0,c3)+A(r1,c2);
oA22 = A(r1,c2)-A(r1,c5)-oA10;
oA23 = A(r2,c5)+oA6+t26;
oA24 = A(r2,c1)+oA3-A(r1,c4);
oA25 = A(r2,c1)-A(r2,c5)-A(r0,c3);
oA26 = A(r1,c1)-oA6-A(r0,c0);
oA27 = oA10-A(r0,c4)-A(r1,c1);
oA28 = A(r2,c1)+t19-t22;
oA30 = t18+t25-A(r1,c1);
oA34 = A(r0,c0)+A(r0,c5)-A(r2,c2)+t27;
oA35 = A(r1,c3)-oA3-A(r0,c1);
oA37 = A(r0,c2)+t27;
oA38 = A(r0,c5)+A(r1,c4)+A(r2,c3);
oA39 = A(r1,c0)+A(r1,c5)-A(r0,c4);
oA1 = A(r0,c4);
oA2 = A(r0,c2);
oA4 = A(r2,c2);
oA5 = A(r0,c3);
oA7 = A(r0,c5);
oA8 = A(r2,c4);
oA11 = A(r1,c4);
oA12 = A(r2,c5);
oA13 = A(r1,c5);
oA14 = A(r0,c1);
oA15 = A(r0,c0);
oA17 = A(r1,c2);
oA19 = A(r2,c0);
oA29 = A(r2,c1);
oA31 = A(r1,c1);
oA32 = A(r2,c3);
oA33 = A(r1,c0);
oA36 = A(r1,c3);

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t18 = B(r1,c0)+B(r3,c2);
t19 = B(r0,c0)-B(r5,c1);
t20 = B(r0,c2)-B(r5,c2);
t21 = B(r2,c1)+B(r4,c0);
t22 = t18+t21;
oB5 = B(r2,c1)-t20;
oB4 = B(r4,c1)+t18;
t25 = B(r3,c1)+B(r0,c1);
oB6 = B(r2,c2)-t19;
oB2 = t25-oB5;
oB10 = B(r2,c0)+t22;
t30 = B(r1,c2)-B(r2,c2);
oB39 = B(r1,c1)+B(r3,c1)-B(r2,c0);
t32 = B(r3,c0)+B(r1,c0);
oB7 = B(r2,c1)+t32;
oB13 = t21-B(r4,c1);
oB17 = B(r2,c1)*2+B(r4,c0)-B(r4,c1)-t20-oB39;
oB18 = B(r4,c2)+B(r5,c1)-B(r2,c1);
oB19 = -B(r3,c0)-B(r4,c2)-oB6;
oB21 = B(r0,c0)+B(r5,c0)-oB2;
oB23 = B(r5,c0)-oB6-B(r3,c1);
oB25 = t25-B(r5,c0)-t19;
oB27 = B(r3,c0)-B(r3,c2)-B(r2,c0);
oB28 = t30-oB4;
oB30 = -B(r0,c2)-B(r3,c0)-B(r3,c1);
oB31 = B(r0,c0)-B(r1,c2)+oB10;
oB32 = B(r0,c1)+B(r1,c0)-B(r5,c2);
oB33 = B(r5,c1)+t22-t30;
oB37 = B(r0,c1)-B(r1,c1)-oB4;
oB38 = B(r4,c2)-t19+t32;
oB0 = B(r2,c2);
oB1 = B(r2,c0);
oB3 = B(r2,c1);
oB8 = B(r5,c1);
oB9 = B(r3,c1);
oB11 = B(r1,c0);
oB12 = B(r0,c1);
oB14 = B(r3,c0);
oB15 = B(r0,c0);
oB16 = B(r1,c1);
oB20 = B(r3,c2);
oB22 = B(r4,c0);
oB24 = B(r4,c1);
oB26 = B(r5,c0);
oB29 = B(r5,c2);
oB34 = B(r1,c2);
oB35 = B(r0,c2);
oB36 = B(r4,c2);

iC0 = FMMa363_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMMa363_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMMa363_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMMa363_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMMa363_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMMa363_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMMa363_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMMa363_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMMa363_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMMa363_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMMa363_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMMa363_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMMa363_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMMa363_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMMa363_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMMa363_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMMa363_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMMa363_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMMa363_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMMa363_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMMa363_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMMa363_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMMa363_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMMa363_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMMa363_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMMa363_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMMa363_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMMa363_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMMa363_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMMa363_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMMa363_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMMa363_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMMa363_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMMa363_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMMa363_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMMa363_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMMa363_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMMa363_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMMa363_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMMa363_mul( oA39, oB39, nmin, peeling, level);
b1 = iC30-iC2;
b2 = iC27+iC37;
v47 = iC20+iC32;
v48 = iC39-iC25;
v52 = iC21+iC35;
b8 = iC0-iC28;
v51 = iC23-iC33;
v55 = iC38-iC24;
b10 = iC37+iC34+iC4-iC10;
b13 = iC18-iC15;
b14 = iC12+iC16;
b15 = iC23+iC6+iC26+iC9;
b16 = iC38-iC32-iC11-iC7;
v46 = iC1+b14;
b26 = iC8-iC9-v48;
b30 = iC4-iC7-v47-v46;
b35 = iC17-iC2+b8+b10;
oC7 = iC36+v55-iC26-b2+v48;
oC6 = iC31-b14-iC15+b10+b15;
oC0 = iC12-b16-b35-b13-iC29;
oC2 = iC22+iC34-v52+v51+v47;
oC1 = iC14+b1+b2+b30;
oC3 = iC13+b8-v51+v46+b26;
oC5 = iC5+iC8-iC21+iC25+b35;
oC8 = iC3+v55+b13+v52-b26+b30;
oC4 = oC6+oC0+iC19+iC0+b1+b16-b15;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
