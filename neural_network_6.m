function [Y_1,Y_2] = neural_network_6(W_1,W_2,W_3,W_4,W_5,W_6,X,multiply)
    Y = multiply(W_1,X);
    Y = ReLu_mat(Y);
    Y = multiply(W_2,Y);
    Y = ReLu_mat(Y);
    Y = multiply(W_3,Y);
    Y = ReLu_mat(Y);
    Y = multiply(W_4,Y);
    Y = ReLu_mat(Y);
    Y = multiply(W_5,Y);
    Y = ReLu_mat(Y);
    Y = multiply(W_6,Y);
    Y = double(Y);
    Y_1 = real(Y);
    Y_2 = imag(Y);
end