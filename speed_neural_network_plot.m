plot(x,T_R,'-k',x,T_G,'-b',x,T_N,'-r');
xlabel('log condition number of weight matrix',extraInputs{:});
ylabel('speed',extraInputs{:});
legend({'Regular', 'Gauss', 'New'},extraInputs{:});
legend('Location','northeast');
title('speed (neural network)','Interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;