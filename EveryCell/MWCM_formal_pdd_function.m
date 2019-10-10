
clc; clear; close;
set(0,'defaultAxesFontSize',24);

p_bar = 7200; % Total protein number
n = 4.615;

% set(gcf, 'Position',  [1200, 400, 800, 640])
set(gcf, 'Position',  [1200, 400, 1040, 832])

p = 0:1:9000;
beta = 0.001*exp(n*p/p_bar) - 0.001;
plot(p,beta)
% title('\beta = 0.001e^{(6.91*p/p_{ss})} - 0.001    for p_{ss}=7200')
xlabel('$p$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
xlim([0 9000])
% yticks([0,0.33,0.66,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0])
grid on
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
