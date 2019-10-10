

set(0,'DefaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',18);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
c1 = [0.85 0.85 0.85];
c2 = [0.55 0.55 0.55];
c3 = [0.1 0.1 0.1];




%% E/ M

set(gcf, 'Position',  [1200, 4, 800, 640])
plot(EM_info(:,1), EM_info(:,2), 'Color', E_yellow)
hold on
plot(EM_info(:,1), EM_info(:,3), 'Color', M_orange)
xlabel('Time / h')
ylabel('Cell count')
legend('E','M','Location','East');
grid on
hold off


%% mRNA
num = 1;

set(gcf, 'Position',  [1200, 600, 800, 640])
plot(pop{num}(:,1), pop{num}(:,2))
hold on
plot(pop{num}(:,1), pop{num}(:,3))
xlabel('Time / h')
ylabel('mRNA count')
l = legend('Genomic','Heterologous','Location','NorthEast');
grid on
hold off

%% Protein

num = 7;

set(gcf, 'Position',  [1800, 400, 1600, 640])
plot(pop{num}(:,1), pop{num}(:,4))
hold on
plot(pop{num}(:,1), pop{num}(:,5))
xlabel('Time / h')
ylabel('Protein count')
ylim([0 7500])
% legend('Genomic','Heterologous','Location','SouthWest');
grid on
hold off


%%

set(gcf, 'Position',  [2100, 400, 800, 640])

histogram(p_mean)
xlim([0 7800])
xlabel('Protein quantity upon division')
ylabel('Frequency')


