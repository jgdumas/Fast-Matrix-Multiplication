x = 1:20;
x = 33 + x;
x = floor(2.^x);
x = log10(x);
extraInputs = {'interpreter','latex','fontsize',14};
plot(x,E_R,'-k',x,E_G,'-b',x,E_N,'-r');
xlabel('log condition number of matrix',extraInputs{:});
ylabel('relative error',extraInputs{:});
legend({'Regular', 'Gauss', 'New'},extraInputs{:});
legend('Location','northeast');
title('accuracy (unitary transform)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;