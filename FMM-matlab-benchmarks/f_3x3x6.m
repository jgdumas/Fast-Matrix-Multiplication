function C = f_3x3x6(A, B, nmin)
%f_3x3x6   Smirnov's 3x3x6 fast matrix multiplication algorithm (Winograd variant).
%          C = f_3x3x6(A, B, NMIN), where A and B are matrices of dimension
%          a power of 2, computes the product C = A*B.
%          This algorithm is used recursively until dimension <= NMIN is reached,
%          at which point standard multiplication is used.
%          The default is NMIN = 3 (which minimizes the total number of
%          operations).

%          Reference:

if nargin < 3, nmin = 3; end

ca = size(A,1);cb = size(B,2);
if ((ca <= nmin) | (cb <=nmin))
   C = A*B;
else
if ((rem(ca,3) ~= 0) | (rem(size(A,2),3) ~=0))
   error('The first matrix dimension must be a multiple of 3.')
end
if ((rem(cb,3) ~= 0) | (rem(size(B,2),6) ~=0))
   error('The second matrix dimension must be a multiple of [3,6].')
end

   atom = ca/3;
   i = 1:atom; j = (atom+1):(2*atom); k=(2*atom+1):ca;
   l = (ca+1):(ca+atom); m=(ca+atom+1):(ca+2*atom); n=(ca+2*atom+1):(cb);

t9=A(i,i)+A(i,k);
t10=A(k,i)+A(k,j);
t11=A(k,i)-A(k,j);
t12=A(j,j)+t9;
t13=A(j,k)-t10;
t14=t12+A(j,k);
t15=A(i,i)+A(i,j);
l_o2=t11+t14;
l_o4=A(j,j)-t9-t13;
l_o12=t14-t11;
l_o13=t12-t13;
l_o16=t9-A(k,i)-A(k,k);
l_o20=t15-t10;
l_o22=t10+t15;
l_o23=A(i,j)+t11-A(i,i);
l_o25=A(j,j)-t11-A(j,i);
v40=l_o4-l_o13;
l_o32=l_o12+l_o23;
v42=l_o22-l_o13;
l_o9=l_o2+v40;
l_o34=l_o20+l_o4;
v45=l_o12-l_o25;
v46=l_o32-l_o2;
v47=l_o23-l_o9;
l_o0=l_o34-l_o23;
v49=l_o16-l_o20;
l_o10=l_o22-l_o34;
v51=v42-v46;
v52=l_o16+l_o9;
l_o7=l_o32-l_o22;
l_o6=l_o23-v42;
l_o3=l_o22-v47;
l_o8=-v51;
l_o21=-v46;
l_o5=l_o12+v40;
l_o15=v46-l_o34;
l_o24=-l_o4-v45;
l_o30=-v52;
l_o36=-l_o25-v51;
l_o39=v40+v45;
l_o1=l_o20-v42;
l_o11=l_o32-l_o20;
l_o29=l_o16-l_o12;
l_o31=-v42-v49;
l_o26=l_o2-v45;
l_o27=l_o10-v45;
l_o14=v47-l_o20;
l_o37=-v45;
l_o35=-v47;
l_o33=-v42;
l_o28=-l_o4-l_o16;
l_o38=l_o0-l_o25;
l_o17=l_o16+v40;
l_o18=v52-l_o7;
l_o19=-l_o23-v49;

t18=B(j,k)+B(k,k);
t20=B(i,m)+B(k,j)/8;
t21=B(j,k)-B(k,k);
t22=B(i,n)+t18;
r23=(B(i,l)-B(j,j))/8;
t23=B(k,i)-r23;
r24=B(k,l)/8;
t24=B(j,m)+r24;
t26=B(i,n)+t21;
r29=B(j,n)/8;
t30=B(i,i)+t18/8;
t31=B(j,i)-r24;
t32=B(k,i)+r23;
t33=B(k,m)-r29+B(i,j)/8;
t34=B(j,m)-B(k,n)/8;
t35=B(k,i)+(B(i,j)-B(j,l))/8;
r36=(B(k,j)-B(k,n))/8;
t36=B(k,m)-r36;
t38=t20+t24;
t39=B(i,m)+t35;
t41=t20-t24;
t43=t20+t30-B(j,j)/8;
r44=B(i,l)/8;
r46=(B(i,k)+B(j,n))/8;
t46=t32+r46;
t47=t23+r46;
t49=t30+r44;
t50=t31+t33;
t51=B(j,i)+B(j,m)-B(j,l)/8-r29;
t52=B(i,i)-r44-t21/8;
t53=B(k,m)+r36;
t54=B(i,m)-t35;
t55=t31-t33;
r56=t22/8;
r57=t26/8;
r58=(B(k,l)-t22)/8;
r59=(B(k,l)+t26)/8;
r_o0=-t41-t46;
r_o1=t38-t47;
r_o2=t41-t46;
r_o3=t38+t47;
r_o8=t55-t52;
r_o9=t50-t49;
r_o10=-t52-t55;
r_o11=t49+t50;
r_o12=t54+r56-t34;
r_o13=r57-t34-t54;
r_o14=-t34-t39-r56;
r_o15=t34-t39+r57;
r_o17=t32+t36+r59;
r_o19=t32+t53-r58;
r_o20=-t43-t51;
r_o23=t43-t51;
r_o29=t36-r58-t23;
r_o31=t53+r59-t23;
v40=r_o0-r_o15;
v41=r_o2-r_o12;
v42=r_o1+r_o13;
r_o6=v40+r_o8;
r_o5=r_o9-v41;
v45=r_o3+r_o14;
v47=r_o23+r_o9;
v48=r_o11-(r_o10-r_o20);
r_o32=r_o11+r_o12-r_o29;
r_o25=r_o12-r_o23-r_o0;
r_o21=r_o0-r_o2+r_o8+v47;
r_o24=-v42-v48;
r_o34=r_o10+r_o15-r_o31;
r_o27=r_o14-r_o20-r_o1;
r_o22=-r_o1-r_o3-v48;
r_o18=-r_o19-v40-v45;
r_o16=v42-r_o17-v41;
r_o30=-r_o29-v41-v45;
r_o28=v42-r_o31-v40;
r_o4=r_o10-v42;
r_o7=r_o11+v45;
r_o33=r_o31+r_o6-r_o1;
r_o26=r_o6+v47;
r_o39=r_o13-r_o17+r_o5;
r_o35=r_o3+r_o29-r_o5;
r_o38=r_o11-r_o19-r_o0;
r_o37=r_o2-r_o10+r_o17;
r_o36=-r_o14-r_o19-r_o6;

prod0=f_3x3x6(l_o0,r_o0,nmin);
prod1=f_3x3x6(l_o1,r_o1,nmin);
prod2=f_3x3x6(l_o2,r_o2,nmin);
prod3=f_3x3x6(l_o3,r_o3,nmin);
prod4=f_3x3x6(l_o4,r_o4,nmin);
prod5=f_3x3x6(l_o5,r_o5,nmin);
prod6=f_3x3x6(l_o6,r_o6,nmin);
prod7=f_3x3x6(l_o7,r_o7,nmin);
prod8=f_3x3x6(l_o8,r_o8,nmin);
prod9=f_3x3x6(l_o9,r_o9,nmin);
prod10=f_3x3x6(l_o10,r_o10,nmin);
prod11=f_3x3x6(l_o11,r_o11,nmin);
prod12=f_3x3x6(l_o12,r_o12,nmin);
prod13=f_3x3x6(l_o13,r_o13,nmin);
prod14=f_3x3x6(l_o14,r_o14,nmin);
prod15=f_3x3x6(l_o15,r_o15,nmin);
prod16=f_3x3x6(l_o16,r_o16,nmin);
prod17=f_3x3x6(l_o17,r_o17,nmin);
prod18=f_3x3x6(l_o18,r_o18,nmin);
prod19=f_3x3x6(l_o19,r_o19,nmin);
prod20=f_3x3x6(l_o20,r_o20,nmin);
prod21=f_3x3x6(l_o21,r_o21,nmin);
prod22=f_3x3x6(l_o22,r_o22,nmin);
prod23=f_3x3x6(l_o23,r_o23,nmin);
prod24=f_3x3x6(l_o24,r_o24,nmin);
prod25=f_3x3x6(l_o25,r_o25,nmin);
prod26=f_3x3x6(l_o26,r_o26,nmin);
prod27=f_3x3x6(l_o27,r_o27,nmin);
prod28=f_3x3x6(l_o28,r_o28,nmin);
prod29=f_3x3x6(l_o29,r_o29,nmin);
prod30=f_3x3x6(l_o30,r_o30,nmin);
prod31=f_3x3x6(l_o31,r_o31,nmin);
prod32=f_3x3x6(l_o32,r_o32,nmin);
prod33=f_3x3x6(l_o33,r_o33,nmin);
prod34=f_3x3x6(l_o34,r_o34,nmin);
prod35=f_3x3x6(l_o35,r_o35,nmin);
prod36=f_3x3x6(l_o36,r_o36,nmin);
prod37=f_3x3x6(l_o37,r_o37,nmin);
prod38=f_3x3x6(l_o38,r_o38,nmin);
prod39=f_3x3x6(l_o39,r_o39,nmin);

b7=prod39+prod35;
b8=prod7-prod32;
v50=prod19+prod23;
b13=prod12-prod0+prod23-prod25;
v52=prod35+prod36;
b15=prod19-prod14-prod6-prod36;
v53=prod35-prod3;
b20=-prod29-prod5-v53;
b25=-prod20-v50;
b27=b25-prod37-prod12-prod4;
v43=-b15;
v46=-b13;
v44=-b20;
b30=b15-b8+v44;
v41=prod29-b27;
z38=prod38+prod32+v46;
z34=prod34+prod37+b13;
z33=prod33+prod36+b7;
z31=prod31+v43-prod27;
z30=prod30+v43-prod26;
z28=prod28+prod25+b27;
z24=prod24+v41;
z18=prod18+b20-prod21;
z17=prod17+prod22-v44;
z16=prod16-b25;
z15=prod15+prod21+v46+b30;
z13=prod13+prod39-prod22+v53+v41;
z11=prod11+prod12-prod32+prod29+v50;
z10=prod10+prod27-prod37-prod22+b8;
z9=prod9+prod6+prod23-prod26+v52;
z8=prod8+prod21-prod25+prod5+v52;
z2=prod2-(prod26+b30);
z1=prod1+prod20-prod27+prod14-b7;
t57=z34+z33;
t53=z34-z33;
t58=z38-z24;
t51=z38+z24;
t63=z28+z18;
b36=z18-z28;
t61=z31-z17;
b37=-z31-z17;
b38=-z30-z16;
t39=z30-z16;
t59=z15-z13;
t44=z15+z13;
b39=z9-z11;
b40=-z11-z9;
t50=z10+z8;
b41=z8-z10;
b42=z1-z2;
b43=-z2-z1;
t60=-b39;
t33=-t59;
t56=-t58;
t25=t60+t58;
t31=t59+b36;
b44=b38-t60;
t29=t63-b37;
t24=t44+t53;
b47=t50-b40+t53;
o7=t51-b42;
b48=t39-t61+b42;
b49=t56-b43;
o5=t33-t57;
t28=t57-b41;
b50=t63+b37+b41;
o15=t61+t56+t39+t31;
b51=t29-b38;
t19=t44+t29;
o3=b39+b42+t28;
o0=(o5-b49)/8;
o8=b41+t33+t25;
o4=(t28+t25)/8;
o6=(t51+b42-t24)/8;
b53=t50+b40+b51;
o2=t53+b43+b51;
o10=(b47-t56)/8;
o17=b44+b49+b50;
o12=(b50-b44)/8;
o14=b36-t53+b48;
o16=(t31-b48)/8;
o9=t51+b38-t19;
o13=b38+b47+t19;
o11=b43+t56+b53;
o1=t24+b53;

   C = [
o0  o1  o2  o3  o4  o5 ;
o6  o7  o8  o9  o10 o11 ;
o12 o13 o14 o15 o16 o17
	   ];

end 
