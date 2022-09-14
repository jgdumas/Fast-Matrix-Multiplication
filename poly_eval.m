function [P_real,P_imag] = poly_eval(X,b,multiply)
    %evaluate polynomial b at X
    P = X;
    n = size(X);
    n = n(1);
    m = length(b)-1;
    S = b(1)*eye(n) + b(2)*X;
    for k = 2:m
        P = multiply(P,X);
        S = S + b(k+1)*P;
    end
    P_real = real(S);
    P_imag = imag(S);
end