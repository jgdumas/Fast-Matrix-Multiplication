x = 1:100;
x = 23 + x;
x = floor(1.24.^x);
x = log10(x);
extraInputs = {'interpreter','latex','fontsize',14};
extraInputs_2 = {'interpreter','latex','fontsize',10};
plot(x(20:100),T_R(20:100),'-k',x(20:100),T_G(20:100),'-b',x(20:100),T_N(20:100),'-r');
xlabel('log condition number of matrix',extraInputs{:});
ylabel('speed (random matrix)',extraInputs{:});
legend({'Regular', 'Gauss', 'New'},extraInputs_2{:});
legend('Location','northeast');
title('speed (matrix multiplication)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;