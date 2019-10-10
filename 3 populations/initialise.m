
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
lllll;

%% Import, make structure

% Master RHS
config_E_CELL;
config_U_CELL;
config_M_CELL;

bioreactor.E_CELL = E_CELL;
bioreactor.U_CELL = U_CELL;
bioreactor.M_CELL = M_CELL;
bioreactor.y_order = {'E_CELL','U_CELL','M_CELL'};

%% Loop

z_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13];
zu_list = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12];
prom_plus_list = [0,1,2,3,4,5,7.5,10,15,20,25,30,40,50,100];
RBS_plus_list = [0,1,2,3,4,5,7.5,10,15,20,25,30,40,50];

list = z_list;

% Store info
GR_delta = zeros(1, length(list));
t_Mall = zeros(1, length(list));
t_Hgone = zeros(1, length(list));

for i = 1:length(list)

%% Parameters

% p.z = 10e-4;
p.z = z_list(i);

% p.zu = 10e-3;
p.zu = zu_list(i);

p.prom_plus = 10;
% p.prom_plus = prom_plus_list(i);

p.RBS_plus = 10;
% p.RBS_plus = RBS_plus_list(i);

p.Ufactor = 0.5;
% p.Ufactor = Ufactor_list(i);

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
Ucell_0 = p.N - Ecell_0;
Mcell_0 = p.N - Ecell_0;

e_0 = 5.46e5; % Energy
TX_0 = [1255, 74, 1343, 74, 74]; % RNApol-DNA complex
m_0 = [301, 1030, 278, 1030, 1030]; % mRNAs
TL_0 = [1988, 400, 1970, 400, 400]; % Ribosome-mRNA complex
p_0 = [9, 25841, 1801, 25841, 25841]; % Proteins

var = [e_0, TX_0(1:4), m_0(1:4), TL_0(1:4), p_0(1:4),...
         TX_0(end), m_0(end), TL_0(end), p_0(end)];     
var_E = [Ecell_0, var];
var_U = [Ucell_0, var];
var_M = [Mcell_0, var];

x0 = [var_E, var_U, var_M];

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

Ucell = Y(:,L*2-21);
e_U = Y(:,L*2-20);
TX_R_U = Y(:,L*2-19);
TX_C_U = Y(:,L*2-18);
TX_P_U = Y(:,L*2-17);
TX_Q_U = Y(:,L*2-16);
m_R_U = Y(:,L*2-15);
m_C_U = Y(:,L*2-14);
m_P_U = Y(:,L*2-13);
m_Q_U = Y(:,L*2-12);
TL_R_U = Y(:,L*2-11);
TL_C_U = Y(:,L*2-10);
TL_P_U = Y(:,L*2-9);
TL_Q_U = Y(:,L*2-8);
R_U = Y(:,L*2-7);
C_U = Y(:,L*2-6);
P_U = Y(:,L*2-5);
Q_U = Y(:,L*2-4);
TX_H_U = Y(:,L*2-3);
m_H_U = Y(:,L*2-2);
TL_H_U = Y(:,L*2-1);
H_U = Y(:,L*2);

Mcell = Y(:,L*3-21);
e_M = Y(:,L*3-20);
TX_R_M = Y(:,L*3-19);
TX_C_M = Y(:,L*3-18);
TX_P_M = Y(:,L*3-17);
TX_Q_M = Y(:,L*3-16);
m_R_M = Y(:,L*3-15);
m_C_M = Y(:,L*3-14);
m_P_M = Y(:,L*3-13);
m_Q_M = Y(:,L*3-12);
TL_R_M = Y(:,L*3-11);
TL_C_M = Y(:,L*3-10);
TL_P_M = Y(:,L*3-9);
TL_Q_M = Y(:,L*3-8);
R_M = Y(:,L*3-7);
C_M = Y(:,L*3-6);
P_M = Y(:,L*3-5);
Q_M = Y(:,L*3-4);
TX_H_M = Y(:,L*3-3);
m_H_M = Y(:,L*3-2);
TL_H_M = Y(:,L*3-1);
H_M = Y(:,L*3);

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

%%%%%%%%%%%%%%%% U_CELL

TX_all_U = TX_R_U + TX_C_U + TX_P_U + TX_Q_U + TX_H_U;
TL_rate_U = p.v_TL .* e_U ./ (p.K_TL + e_U);
TL_all_U = TL_R_U + TL_C_U + TL_P_U + TL_Q_U + TL_H_U;

GR_U = TL_rate_U .* TL_all_U ./ p.mass;

R_proteome_U = R_U .* (p.n_R/3);
C_proteome_U = C_U .* (p.n_C/3);
P_proteome_U = P_U .* (p.n_P/3);
Q_proteome_U = Q_U .* (p.n_Q/3);
H_proteome_U = H_U .* (p.n_H/3);

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

ALLcell = Ecell + Ucell + Mcell;

e_avg = (Ecell.*e_E + Ucell.*e_U + Mcell.*e_M) ./ ALLcell;

m_R_avg = (Ecell.*m_R_E + Ucell.*m_R_U + Mcell.*m_R_M) ./ ALLcell;
m_C_avg = (Ecell.*m_C_E + Ucell.*m_C_U + Mcell.*m_C_M) ./ ALLcell;
m_P_avg = (Ecell.*m_P_E + Ucell.*m_P_U + Mcell.*m_P_M) ./ ALLcell;
m_Q_avg = (Ecell.*m_Q_E + Ucell.*m_Q_U + Mcell.*m_Q_M) ./ ALLcell;
m_H_avg = (Ecell.*m_H_E + Ucell.*m_H_U + Mcell.*m_H_M) ./ ALLcell;

R_avg = (Ecell.*R_E + Ucell.*R_U + Mcell.*R_M) ./ ALLcell;
C_avg = (Ecell.*C_E + Ucell.*C_U + Mcell.*C_M) ./ ALLcell;
P_avg = (Ecell.*P_E + Ucell.*P_U + Mcell.*P_M) ./ ALLcell;
Q_avg = (Ecell.*Q_E + Ucell.*Q_U + Mcell.*Q_M) ./ ALLcell;
H_avg = (Ecell.*H_E + Ucell.*H_U + Mcell.*H_M) ./ ALLcell;

TX_all_avg = (Ecell.*TX_all_E + Ucell.*TX_all_U + Mcell.*TX_all_M) ./ ALLcell;
TL_all_avg = (Ecell.*TL_all_E + Ucell.*TL_all_U + Mcell.*TL_all_M) ./ ALLcell;

GR_avg = (Ecell.*GR_E + Ucell.*GR_U + Mcell.*GR_M) ./ ALLcell;

R_proteome_avg = (Ecell.*R_proteome_E + Ucell.*R_proteome_U + Mcell.*R_proteome_M) ./ ALLcell;
C_proteome_avg = (Ecell.*C_proteome_E + Ucell.*C_proteome_U + Mcell.*C_proteome_M) ./ ALLcell;
P_proteome_avg = (Ecell.*P_proteome_E + Ucell.*P_proteome_U + Mcell.*P_proteome_M) ./ ALLcell;
Q_proteome_avg = (Ecell.*Q_proteome_E + Ucell.*Q_proteome_U + Mcell.*Q_proteome_M) ./ ALLcell;
H_proteome_avg = (Ecell.*H_proteome_E + Ucell.*H_proteome_U + Mcell.*H_proteome_M) ./ ALLcell;



%% Loop calcs

GR_delta(i) = GR_M(end) - GR_E(end);

% Time until all M
for k = 1:length(T)
    if Mcell(k) > p.N - 0.001*p.N
        t_Mall(i) = T(k);
        break
    end    
end

% Time until all H gone
for k = 1:length(T)
    if H_avg(k) < 2
        t_Hgone(i) = T(k);
        break
    end   
end

toc
          
end


           
%% EM plot

close all
set(gcf, 'Position',  [1000, 600, 800, 640])

plot(T, Ecell, 'Color', E_yellow, 'LineWidth', 3)
% plot(T, Ecell, '--', 'Color', E_yellow, 'LineWidth', 3)
hold on
plot(T, Ucell, 'Color', U_orange, 'LineWidth', 3)
% plot(T, Mcell, '--', 'Color', M_red, 'LineWidth', 3)
plot(T, Mcell, 'Color', M_red, 'LineWidth', 3)
% plot(T, Mcell, '--', 'Color', M_red, 'LineWidth', 3)

% legend('E', 'U', 'M')
xlabel('Time / arb. unit')
ylabel('Number of cells')
xlim([0 120])
% xlim([0 180])
ylim([0 p.N])
grid on
axis square



%% Growth Rate plot

% close all
% set(gcf, 'Position',  [1000, 0, 800, 640])
% 
% plot(T, GR_E, 'Color', E_yellow)
% hold on
% plot(T, GR_U, 'Color', U_orange)
% plot(T, GR_M, 'Color', M_red)
% plot(T, GR_avg, 'k--')
% 
% legend('E', 'U', 'M', 'Avg')
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
% plot(T, H_U)
% plot(T, H_M)
% plot(T, H_avg, 'k--')
% legend('E', 'U', 'M', 'Avg')
% xlabel('Time / h')
% ylabel('Relative H')
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

close all
set(gcf, 'Position',  [1000, 600, 800, 640])

% plot(z_list, t_Mall)
% plot(GR_delta, t_Mall, 'LineWidth', 2)
% plot(prom_plus_list, t_Mall)
% plot(RBS_plus_list, t_Mall)

% plot(GR_delta(2:end), t_Hgone(2:end))
% hold on

semilogx(z_list, t_Hgone)
% xlim([1e-13 1])
% ylim([0 350])

% xlabel('\Delta GR between E and M / h^{-1}')
xlabel('Probability of producing M')
ylabel('Time until no H / arb. unit')
% ylabel('Time until all mutant')
% ylabel('Time until no H')
axis square
grid on


