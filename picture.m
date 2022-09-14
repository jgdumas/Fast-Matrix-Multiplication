x = 1:100;
x = 23 + x;
x = floor(1.24.^x);
x = log10(x);
extraInputs = {'interpreter','latex','fontsize',14};
extraInputs_2 = {'interpreter','latex','fontsize',3};
plot(x,E_R,'-k',x,E_G,'-b',x,E_N,'-r');
xlabel('log condition number of matrix',extraInputs{:});
ylabel('accuracy (random matrix)',extraInputs{:});
legend({'Regular', 'Gauss', 'New'},extraInputs_2{:});
legend('Location','northwest');
title('accuracy (matrix multiplication)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;