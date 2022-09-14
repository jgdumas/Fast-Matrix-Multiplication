x = 1:20;
x = 33 + x;
x = floor(2.^x);
x = log10(x);
extraInputs = {'interpreter','latex','fontsize',14};
extraInputs_2 = {'interpreter','latex','fontsize',5};
plot(x,E_R,'-k',x,E_G,'-b',x,E_N,'-r');
xlabel('log condition number of weight matrix',extraInputs{:});
ylabel('relative forward error',extraInputs{:});
legend({'Regular', 'Gauss', 'New'},extraInputs_2{:});
legend('Location','northwest');
title('accuracy (neural network)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;