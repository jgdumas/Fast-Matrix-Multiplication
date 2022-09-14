function [Tr,Tg,Tn] = speed(n)
    Tr = 0;
    Tg = 0;
    Tn = 0;
    for k = 1:10
        A = rand(n,n);
        B = rand(n,n);
        C = rand(n,n);
        D = rand(n,n);
        tic
        [MR_1,MR_2] = f(A,B,C,D);
        Tr = Tr + toc/10;
        tic
        [MG_1,MG_2] = g(A,B,C,D);
        Tg = Tg + toc/10;
        tic
        [MN_1,MN_2] = h(A,B,C,D);
        Tn = Tn + toc/10;
    end
end