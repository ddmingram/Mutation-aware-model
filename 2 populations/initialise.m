
clc; clear;

set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',24);

E_yellow = [237 177 32]/255;
U_orange = [224 135 0]/255;
M_red = [255, 51, 51]/255;
lightgrey = [217 217 217]/255;
mediumgrey = [140 140 140]/255;
darkgrey = [25 25 25]/255;
R_blue = [0 114 189]/255;
G_green = [119 172 48]/255;
H_yellow = [237 177 32]/255;

tic


%% Import, make structure

% Master RHS
config_E_CELL;
config_M_CELL;

bioreactor.E_CELL = E_CELL;
bioreactor.M_CELL = M_CELL;
bioreactor.y_order = {'E_CELL','M_CELL'};

%% Loop

divisions = 10;

z_list = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13];
prom_plus_list = [0,1,2,3,4,5,7.5,10,15,20,25,30,40,50,75,100];
% prom_plus_list = linspace(0, 50, divisions);
RBS_plus_list = linspace(0, 50, divisions);

list = prom_plus_list;

% Store info
% GR_E_store = zeros(divisions); % Depending on prom (j) and RBS (i)
GR_E_store = zeros(1, length(list));
GR_M_store = zeros(1, length(list));
GR_delta = zeros(1, length(list));
t_Mall = zeros(1, length(list));
t_Hgone = zeros(1, length(list));


% for i = 1:length(list)
% for j = 1:length(list)

%% Parameters

p.z = 10e-4;
% p.z = z_list(i);

p.prom_plus = 10;
% p.prom_plus = prom_plus_list(i);

p.RBS_plus = 10;
% p.RBS_plus = RBS_plus_list(i);

% Normal Params
p.mass = 10^8;
p.N = 10^6;
p.s = 10;
p.v_e = 1000;
p.K_e = 100;
p.n_R = 22377;
p.n_C = 900;
p.n_P = 24804;
p.n_Q = 900;
p.kb_TX = 10;
p.ku_TX = 1;
p.v_TX = 216000;
p.K_TX = 500;
p.K_Q = 150000;
p.hQ = 4;
p.m_deg = 12;
p.kb_TL = 10;
p.ku_TL = 1;
p.v_TL = 72000;
p.K_TL = 400;
p.plasmid = 1;
p.n_H = 900;
p.prom_minus = 1;
p.RBS_minus = 1;
p.Mfactor = 0;


%% Initial conditions

Ecell_0 = p.N;
Mcell_0 = p.N - Ecell_0;

e_0 = 5.46e5; % Energy
TX_0 = [1255, 74, 1343, 74, 74]; % RNApol-DNA complex
m_0 = [301, 1030, 278, 1030, 1030]; % mRNAs
TL_0 = [1988, 400, 1970, 400, 400]; % Ribosome-mRNA complex
p_0 = [9, 25841, 1801, 25841, 25841]; % Proteins

var = [e_0, TX_0(1:4), m_0(1:4), TL_0(1:4), p_0(1:4),...
         TX_0(end), m_0(end), TL_0(end), p_0(end)];     
var_E = [Ecell_0, var];
var_M = [Mcell_0, var];

x0 = [var_E, var_M];

%% ODE

tspan = [0,10000];

[T,Y] = ode15s(@(t,y) RHS_master(t, y, p, bioreactor),...
               tspan,...
               x0);

L = length(var_E);
           
Ecell = Y(:,L-21);
e_E = Y(:,L-20);
TX_R_E = Y(:,L-19);
TX_C_E = Y(:,L-18);
TX_P_E = Y(:,L-17);
TX_Q_E = Y(:,L-16);
m_R_E = Y(:,L-15);
m_C_E = Y(:,L-14);
m_P_E = Y(:,L-13);
m_Q_E = Y(:,L-12);
TL_R_E = Y(:,L-11);
TL_C_E = Y(:,L-10);
TL_P_E = Y(:,L-9);
TL_Q_E = Y(:,L-8);
R_E = Y(:,L-7);
C_E = Y(:,L-6);
P_E = Y(:,L-5);
Q_E = Y(:,L-4);
TX_H_E = Y(:,L-3);
m_H_E = Y(:,L-2);
TL_H_E = Y(:,L-1);
H_E = Y(:,L);

Mcell = Y(:,L*2-21);
e_M = Y(:,L*2-20);
TX_R_M = Y(:,L*2-19);
TX_C_M = Y(:,L*2-18);
TX_P_M = Y(:,L*2-17);
TX_Q_M = Y(:,L*2-16);
m_R_M = Y(:,L*2-15);
m_C_M = Y(:,L*2-14);
m_P_M = Y(:,L*2-13);
m_Q_M = Y(:,L*2-12);
TL_R_M = Y(:,L*2-11);
TL_C_M = Y(:,L*2-10);
TL_P_M = Y(:,L*2-9);
TL_Q_M = Y(:,L*2-8);
R_M = Y(:,L*2-7);
C_M = Y(:,L*2-6);
P_M = Y(:,L*2-5);
Q_M = Y(:,L*2-4);
TX_H_M = Y(:,L*2-3);
m_H_M = Y(:,L*2-2);
TL_H_M = Y(:,L*2-1);
H_M = Y(:,L*2);

%% Calculations

%%%%%%%%%%%%%%%% E_CELL

TX_all_E = TX_R_E + TX_C_E + TX_P_E + TX_Q_E + TX_H_E;
TL_rate_E = p.v_TL .* e_E ./ (p.K_TL + e_E);
TL_all_E = TL_R_E + TL_C_E + TL_P_E + TL_Q_E + TL_H_E;

GR_E = TL_rate_E .* TL_all_E ./ p.mass;

R_proteome_E = R_E .* (p.n_R/3);
C_proteome_E = C_E .* (p.n_C/3);
P_proteome_E = P_E .* (p.n_P/3);
Q_proteome_E = Q_E .* (p.n_Q/3);
H_proteome_E = H_E .* (p.n_H/3);

%%%%%%%%%%%%%%%% M_CELL

TX_all_M = TX_R_M + TX_C_M + TX_P_M + TX_Q_M + TX_H_M;
TL_rate_M = p.v_TL .* e_M ./ (p.K_TL + e_M);
TL_all_M = TL_R_M + TL_C_M + TL_P_M + TL_Q_M + TL_H_M;

GR_M = TL_rate_M .* TL_all_M ./ p.mass;

R_proteome_M = R_M .* (p.n_R/3);
C_proteome_M = C_M .* (p.n_C/3);
P_proteome_M = P_M .* (p.n_P/3);
Q_proteome_M = Q_M .* (p.n_Q/3);
H_proteome_M = H_M .* (p.n_H/3);

%%%%%%%%%%%%%%%% AVERAGING

ALLcell = Ecell + Mcell;

e_avg = (Ecell.*e_E + Mcell.*e_M) ./ ALLcell;

m_R_avg = (Ecell.*m_R_E + Mcell.*m_R_M) ./ ALLcell;
m_C_avg = (Ecell.*m_C_E + Mcell.*m_C_M) ./ ALLcell;
m_P_avg = (Ecell.*m_P_E + Mcell.*m_P_M) ./ ALLcell;
m_Q_avg = (Ecell.*m_Q_E + Mcell.*m_Q_M) ./ ALLcell;
m_H_avg = (Ecell.*m_H_E + Mcell.*m_H_M) ./ ALLcell;

R_avg = (Ecell.*R_E + Mcell.*R_M) ./ ALLcell;
C_avg = (Ecell.*C_E + Mcell.*C_M) ./ ALLcell;
P_avg = (Ecell.*P_E + Mcell.*P_M) ./ ALLcell;
Q_avg = (Ecell.*Q_E + Mcell.*Q_M) ./ ALLcell;
H_avg = (Ecell.*H_E + Mcell.*H_M) ./ ALLcell;

TX_all_avg = (Ecell.*TX_all_E + Mcell.*TX_all_M) ./ ALLcell;
TL_all_avg = (Ecell.*TL_all_E + Mcell.*TL_all_M) ./ ALLcell;

GR_avg = (Ecell.*GR_E + Mcell.*GR_M) ./ ALLcell;

R_proteome_avg = (Ecell.*R_proteome_E + Mcell.*R_proteome_M) ./ ALLcell;
C_proteome_avg = (Ecell.*C_proteome_E + Mcell.*C_proteome_M) ./ ALLcell;
P_proteome_avg = (Ecell.*P_proteome_E + Mcell.*P_proteome_M) ./ ALLcell;
Q_proteome_avg = (Ecell.*Q_proteome_E + Mcell.*Q_proteome_M) ./ ALLcell;
H_proteome_avg = (Ecell.*H_proteome_E + Mcell.*H_proteome_M) ./ ALLcell;



%% Loop calcs

% GR_delta(i) = GR_M(end) - GR_E(end);
% 
% % Time until all mutated
% for k = 1:length(T)
%     if Ecell(k) < 0.001*p.N
%         t_Mall(i) = T(k);
%         break
%     end    
% end
% 
% % Time until all H gone
% for k = 1:length(T)
%     if H_avg(k) < 2
%         t_Hgone(i) = T(k);
%         break
%     end   
% end
% 
% % GR_E_store(i,j) = GR_E(end);
% GR_E_store(i) = GR_E(end);
% GR_M_store(i) = GR_M(end);

toc

% end
% end
           
     

%% EM plot

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
% 
% plot(T, Ecell, 'Color', E_yellow, 'LineWidth', 3)
% % plot(T, Ecell, '--', 'Color', E_yellow, 'LineWidth', 3)
% hold on
% plot(T, Mcell, 'Color', M_red, 'LineWidth', 3)
% % plot(T, Mcell, '--', 'Color', M_red, 'LineWidth', 3)
% 
% % legend('E', 'M')
% xlabel('Time / arb. unit')
% ylabel('Number of cells')
% % xlim([0 60])
% xlim([0 180])
% ylim([0 1e6])
% grid on
% axis square



%% Growth Rate plot

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
% 
% plot(T, GR_E, 'Color', E_yellow)
% hold on
% plot(T, GR_M, 'Color', M_red)
% plot(T, GR_avg, 'k--')
% 
% legend('E', 'M', 'Average')
% xlabel('Time / h')
% ylabel('Growth Rate')
% xlim([0 50])
% grid on
% axis square



%% H plot

% close all
% set(gcf, 'Position',  [1000, 0, 800, 640])
% 
% plot(T, H_E)
% hold on
% plot(T, H_M)
% legend('E', 'M')
% xlabel('Time / h')
% ylabel('Number of cells')
% xlim([0 50])
% grid on
% axis square



%% mRNA

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
% 
% plot(T, m_R_avg, 'Color', mediumgrey)
% hold on
% plot(T, m_C_avg, 'Color', mediumgrey)
% plot(T, m_P_avg, 'Color', mediumgrey)
% plot(T, m_Q_avg, 'Color', mediumgrey)
% plot(T, m_H_avg, 'Color', H_yellow, 'LineWidth', 3)
% % legend('E', 'M')
% xlabel('Time / arb. unit')
% ylabel('Relative mRNA')
% xlim([0 60])
% axis square
% grid on



%% Proteome

% close all
% set(gcf, 'Position',  [2000, 600, 800, 640])
% 
% all_proteome_avg = R_proteome_avg + C_proteome_avg + P_proteome_avg + Q_proteome_avg + H_proteome_avg;
% TX_proteome_avg = TX_all_avg.*(p.n_P/3);
% TL_proteome_avg = TL_all_avg.*(p.n_R/3);
% 
% plot(T, (R_proteome_avg + TL_proteome_avg)./p.mass, 'Color', mediumgrey)
% hold on
% plot(T, (C_proteome_avg)./p.mass, 'Color', mediumgrey)
% plot(T, (P_proteome_avg + TX_proteome_avg)./p.mass, 'Color', mediumgrey)
% plot(T, (Q_proteome_avg)./p.mass, 'Color', mediumgrey)
% plot(T, (H_proteome_avg)./p.mass, 'Color', H_yellow, 'LineWidth', 3)
% 
% % plot(T, all_proteome_avg + TX_proteome_avg + TL_proteome_avg, 'k--')
% 
% xlabel('Time / arb. unit')
% ylabel('Fraction of proteome')
% xlim([0 60])
% axis square
% grid on



%% Heatmaps

% close all
% 
% heat_GR = flipud(GR_E_store); % Inverts rows
% xvalues = prom_plus_list;
% yvalues = fliplr(RBS_plus_list); % Inverts vector
% 
% h1 = heatmap(xvalues, yvalues, heat_GR);
% h1.xlabel('Promoter Strength')
% h1.ylabel('RBS strength')
% h1.GridVisible = 'off';
% 
% colormap default
% h1.FontSize = 24;



%% Loop plotting

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
% 
% % semilogx(z_list, GR_E_store, 'o', 'Color', E_yellow)
% % hold on
% % semilogx(z_list, GR_M_store, 'o', 'Color', M_red)
% % xlim([1e-13 1])
% % ylim([3.6 4])
% % set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
% 
% plot(GR_delta(2:end), t_Hgone(2:end))
% xlim([0 1])
% xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
% set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
% 
% % semilogx(z_list, t_Hgone, 'Color', [0, 0.4470, 0.7410])
% % xlim([1e-13 1])
% % ylim([0 350])
% 
% xlabel('\Delta GR between E and M / h^{-1}')
% % xlabel('Probability of producing M')
% 
% % ylabel('Growth rate / h^{-1}')
% % ylabel('Time until x / arb. unit')
% % ylabel('Time until all mutant')
% ylabel('Time until no H / arb. unit')
% 
% axis square
% grid on



plot(T, (H_proteome_avg)./p.mass, '--', 'Color', H_yellow)


