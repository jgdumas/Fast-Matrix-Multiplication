function y = relu(x)
    y_real = max(real(x),0);
    y_imag = max(imag(x),0);
    y = y_real + 1j*y_imag;
end