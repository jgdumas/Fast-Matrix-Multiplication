function T = gen_mat_SW_axb(n,b)
    T = 2*rand(n)-1;
    %T = randn(n); use this line if normal distribution is used
    q = idivide(int32(n),int32(b));
    r = rem(n,b);
    if r(1)~=0 q(1)=q(1)+1; end;
    if r(2)~=0 q(2)=q(2)+1; end;
    T = [T,zeros(n(1),b(2)*q(2)-n(2));zeros(b(1)*q(1)-n(1),b(2)*q(2))];
end
