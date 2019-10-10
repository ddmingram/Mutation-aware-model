
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

%% Import

% Master RHS
config_E_CELL;
config_MP_CELL;
config_MR_CELL;
config_MPR_CELL;

bioreactor.E_CELL = E_CELL;
bioreactor.MP_CELL = MP_CELL;
bioreactor.MR_CELL = MR_CELL;
bioreactor.MPR_CELL = MPR_CELL;
bioreactor.y_order = {'E_CELL','MP_CELL','MR_CELL','MPR_CELL'};

%% Loop

zp_list = [1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.001, 0.0001];
prom_plus_list = [0,1,2,3,4,5,7.5,10,15,20,25,30,40,50];
RBS_plus_list = [0,1,2,3,4,5,7.5,10,15,20,25,30,40,50];

list = zp_list;

GR_delta = zeros(1, length(list));
t_finish = zeros(1, length(list));

% for i = 1:length(list)

%% Parameters

p.zp = 5e-4;
% p.zp = zp_list(i);

p.zr = 5e-4;
% p.zr = zr_list(i);

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
Mpcell_0 = p.N - Ecell_0;
Mrcell_0 = p.N - Ecell_0;
Mprcell_0 = p.N - Ecell_0;

e_0 = 5.46e5; % Energy
TX_0 = [1255, 74, 1343, 74, 74]; % RNApol-DNA complex
m_0 = [301, 1030, 278, 1030, 1030]; % mRNAs
TL_0 = [1988, 400, 1970, 400, 400]; % Ribosome-mRNA complex
p_0 = [9, 25841, 1801, 25841, 25841]; % Proteins

var = [e_0, TX_0(1:4), m_0(1:4), TL_0(1:4), p_0(1:4),...
         TX_0(end), m_0(end), TL_0(end), p_0(end)];     
var_E = [Ecell_0, var];
var_Mp = [Mpcell_0, var];
var_Mr = [Mrcell_0, var];
var_Mpr = [Mprcell_0, var];

x0 = [var_E, var_Mp, var_Mr, var_Mpr];

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

Mpcell = Y(:,L*2-21);
e_Mp = Y(:,L*2-20);
TX_R_Mp = Y(:,L*2-19);
TX_C_Mp = Y(:,L*2-18);
TX_P_Mp = Y(:,L*2-17);
TX_Q_Mp = Y(:,L*2-16);
m_R_Mp = Y(:,L*2-15);
m_C_Mp = Y(:,L*2-14);
m_P_Mp = Y(:,L*2-13);
m_Q_Mp = Y(:,L*2-12);
TL_R_Mp = Y(:,L*2-11);
TL_C_Mp = Y(:,L*2-10);
TL_P_Mp = Y(:,L*2-9);
TL_Q_Mp = Y(:,L*2-8);
R_Mp = Y(:,L*2-7);
C_Mp = Y(:,L*2-6);
P_Mp = Y(:,L*2-5);
Q_Mp = Y(:,L*2-4);
TX_H_Mp = Y(:,L*2-3);
m_H_Mp = Y(:,L*2-2);
TL_H_Mp = Y(:,L*2-1);
H_Mp = Y(:,L*2);

Mrcell = Y(:,L*3-21);
e_Mr = Y(:,L*3-20);
TX_R_Mr = Y(:,L*3-19);
TX_C_Mr = Y(:,L*3-18);
TX_P_Mr = Y(:,L*3-17);
TX_Q_Mr = Y(:,L*3-16);
m_R_Mr = Y(:,L*3-15);
m_C_Mr = Y(:,L*3-14);
m_P_Mr = Y(:,L*3-13);
m_Q_Mr = Y(:,L*3-12);
TL_R_Mr = Y(:,L*3-11);
TL_C_Mr = Y(:,L*3-10);
TL_P_Mr = Y(:,L*3-9);
TL_Q_Mr = Y(:,L*3-8);
R_Mr = Y(:,L*3-7);
C_Mr = Y(:,L*3-6);
P_Mr = Y(:,L*3-5);
Q_Mr = Y(:,L*3-4);
TX_H_Mr = Y(:,L*3-3);
m_H_Mr = Y(:,L*3-2);
TL_H_Mr = Y(:,L*3-1);
H_Mr = Y(:,L*3);

Mprcell = Y(:,L*4-21);
e_Mpr = Y(:,L*4-20);
TX_R_Mpr = Y(:,L*4-19);
TX_C_Mpr = Y(:,L*4-18);
TX_P_Mpr = Y(:,L*4-17);
TX_Q_Mpr = Y(:,L*4-16);
m_R_Mpr = Y(:,L*4-15);
m_C_Mpr = Y(:,L*4-14);
m_P_Mpr = Y(:,L*4-13);
m_Q_Mpr = Y(:,L*4-12);
TL_R_Mpr = Y(:,L*4-11);
TL_C_Mpr = Y(:,L*4-10);
TL_P_Mpr = Y(:,L*4-9);
TL_Q_Mpr = Y(:,L*4-8);
R_Mpr = Y(:,L*4-7);
C_Mpr = Y(:,L*4-6);
P_Mpr = Y(:,L*4-5);
Q_Mpr = Y(:,L*4-4);
TX_H_Mpr = Y(:,L*4-3);
m_H_Mpr = Y(:,L*4-2);
TL_H_Mpr = Y(:,L*4-1);
H_Mpr = Y(:,L*4);

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

%%%%%%%%%%%%%%%% Mp_CELL

TX_all_Mp = TX_R_Mp + TX_C_Mp + TX_P_Mp + TX_Q_Mp + TX_H_Mp;
TL_rate_Mp = p.v_TL .* e_Mp ./ (p.K_TL + e_Mp);
TL_all_Mp = TL_R_Mp + TL_C_Mp + TL_P_Mp + TL_Q_Mp + TL_H_Mp;

GR_Mp = TL_rate_Mp .* TL_all_Mp ./ p.mass;

R_proteome_Mp = R_Mp .* (p.n_R/3);
C_proteome_Mp = C_Mp .* (p.n_C/3);
P_proteome_Mp = P_Mp .* (p.n_P/3);
Q_proteome_Mp = Q_Mp .* (p.n_Q/3);
H_proteome_Mp = H_Mp .* (p.n_H/3);

%%%%%%%%%%%%%%%% Mr_CELL

TX_all_Mr = TX_R_Mr + TX_C_Mr + TX_P_Mr + TX_Q_Mr + TX_H_Mr;
TL_rate_Mr = p.v_TL .* e_Mr ./ (p.K_TL + e_Mr);
TL_all_Mr = TL_R_Mr + TL_C_Mr + TL_P_Mr + TL_Q_Mr + TL_H_Mr;

GR_Mr = TL_rate_Mr .* TL_all_Mr ./ p.mass;

R_proteome_Mr = R_Mr .* (p.n_R/3);
C_proteome_Mr = C_Mr .* (p.n_C/3);
P_proteome_Mr = P_Mr .* (p.n_P/3);
Q_proteome_Mr = Q_Mr .* (p.n_Q/3);
H_proteome_Mr = H_Mr .* (p.n_H/3);

%%%%%%%%%%%%%%%% Mpr_CELL

TX_all_Mpr = TX_R_Mpr + TX_C_Mpr + TX_P_Mpr + TX_Q_Mpr + TX_H_Mpr;
TL_rate_Mpr = p.v_TL .* e_Mpr ./ (p.K_TL + e_Mpr);
TL_all_Mpr = TL_R_Mpr + TL_C_Mpr + TL_P_Mpr + TL_Q_Mpr + TL_H_Mpr;

GR_Mpr = TL_rate_Mpr .* TL_all_Mpr ./ p.mass;

R_proteome_Mpr = R_Mpr .* (p.n_R/3);
C_proteome_Mpr = C_Mpr .* (p.n_C/3);
P_proteome_Mpr = P_Mpr .* (p.n_P/3);
Q_proteome_Mpr = Q_Mpr .* (p.n_Q/3);
H_proteome_Mpr = H_Mpr .* (p.n_H/3);

%%%%%%%%%%%%%%%% AVERAGING

ALLcell = Ecell + Mpcell + Mrcell + Mprcell;

e_avg = (Ecell.*e_E + Mpcell.*e_Mp + Mrcell.*e_Mr + Mprcell.*e_Mpr) ./ ALLcell;

m_R_avg = (Ecell.*m_R_E + Mpcell.*m_R_Mp + Mrcell.*m_R_Mr + Mprcell.*m_R_Mpr) ./ ALLcell;
m_C_avg = (Ecell.*m_C_E + Mpcell.*m_C_Mp + Mrcell.*m_C_Mr + Mprcell.*m_C_Mpr) ./ ALLcell;
m_P_avg = (Ecell.*m_P_E + Mpcell.*m_P_Mp + Mrcell.*m_P_Mr + Mprcell.*m_P_Mpr) ./ ALLcell;
m_Q_avg = (Ecell.*m_Q_E + Mpcell.*m_Q_Mp + Mrcell.*m_Q_Mr + Mprcell.*m_Q_Mpr) ./ ALLcell;
m_H_avg = (Ecell.*m_H_E + Mpcell.*m_H_Mp + Mrcell.*m_H_Mr + Mprcell.*m_H_Mpr) ./ ALLcell;

R_avg = (Ecell.*R_E + Mpcell.*R_Mp + Mrcell.*R_Mr + Mprcell.*R_Mpr) ./ ALLcell;
E_avg = (Ecell.*C_E + Mpcell.*C_Mp + Mrcell.*C_Mr + Mprcell.*C_Mpr) ./ ALLcell;
P_avg = (Ecell.*P_E + Mpcell.*P_Mp + Mrcell.*P_Mr + Mprcell.*P_Mpr) ./ ALLcell;
Q_avg = (Ecell.*Q_E + Mpcell.*Q_Mp + Mrcell.*Q_Mr + Mprcell.*Q_Mpr) ./ ALLcell;
H_avg = (Ecell.*H_E + Mpcell.*H_Mp + Mrcell.*H_Mr + Mprcell.*H_Mpr) ./ ALLcell;

TX_all_avg = (Ecell.*TX_all_E + Mpcell.*TX_all_Mp + Mrcell.*TX_all_Mr + Mprcell.*TX_all_Mpr) ./ ALLcell;
TL_all_avg = (Ecell.*TL_all_E + Mpcell.*TL_all_Mp + Mrcell.*TL_all_Mr + Mprcell.*TL_all_Mpr) ./ ALLcell;

GR_avg = (Ecell.*GR_E + Mpcell.*GR_Mp + Mrcell.*GR_Mr + Mprcell.*GR_Mpr) ./ ALLcell;

R_proteome_avg = (Ecell.*R_proteome_E + Mpcell.*R_proteome_Mp + Mrcell.*R_proteome_Mr + Mprcell.*R_proteome_Mpr) ./ ALLcell;
C_proteome_avg = (Ecell.*C_proteome_E + Mpcell.*C_proteome_Mp + Mrcell.*C_proteome_Mr + Mprcell.*C_proteome_Mpr) ./ ALLcell;
P_proteome_avg = (Ecell.*P_proteome_E + Mpcell.*P_proteome_Mp + Mrcell.*P_proteome_Mr + Mprcell.*P_proteome_Mpr) ./ ALLcell;
Q_proteome_avg = (Ecell.*Q_proteome_E + Mpcell.*Q_proteome_Mp + Mrcell.*Q_proteome_Mr + Mprcell.*Q_proteome_Mpr) ./ ALLcell;
H_proteome_avg = (Ecell.*H_proteome_E + Mpcell.*H_proteome_Mp + Mrcell.*H_proteome_Mr + Mprcell.*H_proteome_Mpr) ./ ALLcell;



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

toc

% end


           
%% EM plot

close all
set(gcf, 'Position',  [1000, 600, 800, 640])

plot(T, Ecell, 'Color', E_yellow)
% plot(T, Ecell, '--', 'Color', E_yellow)
hold on
plot(T, Mpcell, 'Color', mediumgrey)
% plot(T, Mpcell, '--', 'Color', mediumgrey)
plot(T, Mrcell, 'Color', darkgrey)
% plot(T, Mrcell, '--', 'Color', darkgrey)
plot(T, Mprcell, 'Color', M_red)
% plot(T, Mprcell, '--', 'Color', M_red)

% legend('E','Mp','Mr', 'Mpr')
xlabel('Time / arb. unit')
ylabel('Number of cells')
xlim([0 1200])
ylim([0 p.N])
grid on
axis square



%% Growth Rate plot

% close all
% set(gcf, 'Position',  [1000, 0, 800, 640])
% 
% plot(T, GR_E, 'Color', E_yellow)
% hold on
% plot(T, GR_Mp, 'Color', mediumgrey)
% plot(T, GR_Mr, 'Color', darkgrey)
% plot(T, GR_Mpr, 'Color', M_red)
% plot(T, GR_avg, 'k--')
% 
% legend('E','Mp','Mr', 'Mpr', 'Average')
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
% plot(T, H_Mp)
% plot(T, H_Mr)
% plot(T, H_Mp)

% legend('E','Mp','Mr', 'Mpr')
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
% 
% xlabel('Time / arb. unit')
% ylabel('Relative mRNA')
% xlim([0 60])
% axis square
% grid on



%% Proteome

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
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
% xlabel('Time / arb. unit')
% ylabel('Fraction of proteome')
% xlim([0 120])
% axis square
% grid on



%% Loop plotting

% close all
% set(gcf, 'Position',  [1000, 600, 800, 640])
% 
% % plot(z_list, t_Mall)
% plot(GR_delta, t_Mall, 'LineWidth', 2)
% % plot(prom_plus_list, t_Mall)
% % plot(RBS_plus_list, t_Mall)
% 
% hold on
% plot(GR_delta(2:end), t_Hgone(2:end), 'LineWidth', 2)
% 
% xlabel('\Delta GR between E and M')
% % xlabel('Mutation probability')
% ylabel('Time until x')
% % ylabel('Time until all mutant')
% % ylabel('Time until no H')
% axis square
% grid on


