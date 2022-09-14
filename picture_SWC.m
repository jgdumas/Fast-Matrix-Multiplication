x = 1:t;
x = 20 + 2*x;
x = log10(x);
extraInputs = {'interpreter','latex','fontsize',14};
%extraInputs_2 = {'interpreter','latex','fontsize',13};
plot(x,E_S,'-r',x,E_W,'-b',x,E_R,'-g');
xlabel('Logarithmic Size of Dimension',extraInputs{:});
ylabel('Relative Error',extraInputs{:});
legend({'Strassen', 'Winograd','Conventional'},extraInputs{:});
legend('Location','northwest');
title('accuracy (uniform random complex matrix)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;