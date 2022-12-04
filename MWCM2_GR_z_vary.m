
%% Set-up

tic % Start timer
close; clc; clear;
tspan = [0,1e5];

set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',28);

                          
%% Collate parameters

a = 1;                  % Change to 60 to convert to /min

n_Int = 0;              % No. intermediate cells
n = n_Int+2;
N = 1e9;                % Value will be updated in core script
mass = 10^8;
nut = 1e4;              % Extracellular nutrient
nq = 1;                 % Nutrient Quality
v_e = 38700/a;          % /h ...= 1/((1/trans)+(1/met))
K_e = 500;              % molecs/cell
n_R = 7459;
n_C = 300;
n_Q = 300;
n_H = 300;
v_TX_R = 55800/a;       % molecs/h/cell
v_TX_C = 248.4/a;         % molecs/h/cell
v_TX_Q = 56940/a;       % molecs/h/cell
K_TX_R = 427;           % molecs/cell
K_TX_nR = 4.38;            % molecs/cell
K_Q = 152000;           % molecs/cell
h_Q = 4;
m_deg = 60*(log(2)/2)/a;       % /h
kb_TL = 60/a;           % cell/h/molecs
ku_TL = 60/a;           % /h
v_TL = 72000/a;         % aa/h (=20 aa/s)
K_TL = 7;               % molecs/cell (=v_TL/K_P = 72000/(180*60))


%% Initialise and store variables

subpop_0 = [N, zeros(1,n-1)]; % Start with all first population (E)

%%% Assign SS values - if H starts at 0
e_0 = ones(1,n)*13.3; 
m_R_0 = ones(1,n)*39.7;
m_C_0 = ones(1,n)*6.6;
m_Q_0 = ones(1,n)*365.7;
m_H_0 = ones(1,n)*0;
TL_R_0 = ones(1,n)*604;
TL_C_0 = ones(1,n)*30.9;
TL_Q_0 = ones(1,n)*1719.6;
TL_H_0 = ones(1,n)*0;
R_0 = ones(1,n)*17.1;
C_0 = ones(1,n)*3633.2;
Q_0 = ones(1,n)*20205;
H_0 = ones(1,n)*0;

var = [subpop_0, e_0,...
       m_R_0, m_C_0, m_Q_0, m_H_0,...
       TL_R_0, TL_C_0, TL_Q_0, TL_H_0,...
       R_0, C_0, Q_0, H_0];

%%% Allocate variables to store results

trials = 2e3;
% To generate GRratio(1->0.2, spaced by 0.025)
% prom_range = [0	136.1	306.2	465.0	646.3	825.0	1020.5	1230.0	1462.7	1700.9	1973.0	2262.1	2568.3	2925.5	3299.6	3724.9	4184.1	4694.3	5272.6	5910.0	6633.3	7432.7	8360.0	9405.7	10613.3	12008.0	13640.8	15562.8	17858.9	20631.3	24033.0	28302.2	33812.9];

divisions = 10; % 11 for heatmaps, 21 for contour maps
prom_range = logspace(2.5, 5, divisions);
z_range = logspace(-12, -3, divisions);

tmut_range = [0.1, linspace(1,48,48), tspan(end)];
Hdrop_range = [0.95, 0.90, 0.75, 0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00001];

t_Hdrop_store = zeros(length(z_range), length(prom_range), length(Hdrop_range));
Harea_Hdrop_store = zeros(length(z_range), length(prom_range), length(Hdrop_range));
Harea_tmut_store = zeros(length(z_range), length(prom_range), length(tmut_range));

track = 1;

for COL = 1:length(z_range)
% for COL = 1:1

z_M = z_range(COL);
% z_M = 1e-9;

for ROW = 1:length(prom_range)
% for ROW = 1:1

prom_E = prom_range(ROW);
% prom_E = 0;

% for k = 1:trials
% prom_E = prom_range(k);

% for k = 1:1
% prom_E = 1e2;

%%% Test
% z_M = z_range(1);
% prom_E = prom_range(1);

prom_I = prom_E/10;
prom_M = 0;
prom_vec = [prom_E, ones(1,n_Int)*prom_I, prom_M];

z_I = 1e-3;
z_vec = [ones(1,n_Int)*z_I, z_M];

params = [n, N, mass,...
          nut, nq, v_e, K_e, n_R, n_C, n_Q, n_H,...
          v_TX_R, v_TX_C, v_TX_Q, K_TX_R, K_TX_nR, K_Q, h_Q, m_deg,...
          kb_TL, ku_TL, v_TL, K_TL,...
          prom_vec, z_vec];


%% 5. Call ODE solver

Opt_1 = odeset('RelTol',1e-6,'AbsTol',1e-9);

[t,y] = ode15s( @(t,y) MWCM2_odes_GR_z_vary(t, y, params),...
               tspan,...
               var,...
               Opt_1);


           
%% 6. Extract variables

x=1;

subpop = y(:, x:x+n-1); x=x+n;
e = y(:, x:x+n-1); x=x+n;
m_R = y(:, x:x+n-1); x=x+n;
m_C = y(:, x:x+n-1); x=x+n;
m_Q = y(:, x:x+n-1); x=x+n;
m_H = y(:, x:x+n-1); x=x+n;
TL_R = y(:, x:x+n-1); x=x+n;
TL_C = y(:, x:x+n-1); x=x+n;
TL_Q = y(:, x:x+n-1); x=x+n;
TL_H = y(:, x:x+n-1); x=x+n;
R = y(:, x:x+n-1); x=x+n;
C = y(:, x:x+n-1); x=x+n;
Q = y(:, x:x+n-1); x=x+n;
H = y(:, x:x+n-1);

% TL_all = TL_R + TL_C + TL_Q + TL_H;
% TL_rate = (v_TL * e) ./ (K_TL + e);
% GR = TL_rate .* TL_all / mass;
% GR_ratio_store(k) = GR(end,1) / GR(end,end);
% disp(k)
% prom_range = prom_range';
% GR_ratio_store = GR_ratio_store';
% end


%% 7. Calculations

% Building blocks
epsilon = nq*C*v_e*nut ./ (K_e + nut);

% Transcription
IQ = (1 ./ (1 + (Q/K_Q).^h_Q));         % AutoInhibition of m_Q

w_R = (v_TX_R * e) ./ (K_TX_R + e);
w_C = (v_TX_C * e) ./ (K_TX_R + e);
w_Q = (v_TX_Q * e) ./ (K_TX_R + e) .* IQ;
w_H = (v_TX_R * e) ./ (K_TX_R + e);

% Translation
TL_rate = (v_TL * e) ./ (K_TL + e);
gamma_R = TL_R .* TL_rate / n_R;
gamma_C = TL_C .* TL_rate / n_C;
gamma_Q = TL_Q .* TL_rate / n_Q;
gamma_H = TL_H .* TL_rate / n_H;

% Summations
m_all = m_R + m_C + m_Q + m_H;
TL_all = TL_R + TL_C + TL_Q + TL_H;
gamma_all = gamma_R + gamma_C + gamma_Q + gamma_H;

% Growth and Dilution
GR = TL_rate .* TL_all / mass;
buffer = sum(subpop,2) - N;

%%%%%%%%%%%%%%%% AVERAGES

e_avg = sum(e.*subpop, 2) ./ sum(subpop,2);
R_avg = sum(R.*subpop, 2) ./ sum(subpop,2);
C_avg = sum(C.*subpop, 2) ./ sum(subpop,2);
Q_avg = sum(Q.*subpop, 2) ./ sum(subpop,2);
TL_all_avg = sum(TL_all.*subpop, 2) ./ sum(subpop,2);
GR_avg = sum(GR.*subpop, 2) ./ sum(subpop,2);

Hrate_pop = sum(GR .* H .* subpop, 2);
Hrate_percell = Hrate_pop ./ N;

Hyield_cumu_pop = cumtrapz(t, Hrate_pop ./ log(2));
Hyield_cumu_percell = Hyield_cumu_pop ./ N;

f_H_pop = Hrate_pop ./ log(2);
f_H_percell = f_H_pop ./ N;

Hyield_pop = trapz(t, Hrate_pop ./ log(2));
Hyield_percell = Hyield_pop ./ N;


%% Test

% close
% plot(t, sum(H.*subpop, 2))
% xlim([0 400])
% 
% hold on

% close
% set(gcf, 'Position',  [600, 200, 650, 500])
% set(0,'defaultAxesFontSize',28);
% set(0,'DefaultLineLineWidth',4);
% 
% plot(t, Hrate_percell ./ log(2), 'k')
% xlim([0 70])
% xticks([0 35 70])
% yticks([0 5e4 10e4])
% xlabel('Time/h')
% ylabel('per cell/h^{−1}')
% axis square

% plot(t, Hrate_percell)
% xlim([0 35])
% xticks([0 15 30])
% yticks([0 5e4 1e5])
% xlabel('Time/h')
% ylabel('H_{rate} per cell/h^{−1}')
% yyaxis right
% plot(t, Hyield_cumu_percell)
% ylabel('Cumulative H_{yield} per cell')
% yyaxis left
% axis square


%% 

% Hrate(1,:) = 0; % t0 produces NaN, hence need to convert to 0
% Hyield_percell = sum(Hrate.*subpop, 2) ./ sum(subpop,2);

%%%%%%%%%%%%%%%% Various measurement metrics

[Hmax, Hmax_t_idx] = max(f_H_percell);

% Time after H drops to % of max
for zstack = 1:length(Hdrop_range)
    for i = 1:length(f_H_percell)
        if (i > Hmax_t_idx) && (f_H_percell(i) <= Hdrop_range(zstack)*Hmax)
            t_Hdrop_store(ROW,COL,zstack) = t(i);
            break
        end
    end
end

% Hyield after H drops to % of max
for zstack = 1:length(Hdrop_range)
    for i = 1:length(f_H_percell)
        if (i > Hmax_t_idx) && (f_H_percell(i) <= Hdrop_range(zstack)*Hmax)
            Harea_Hdrop_store(ROW,COL,zstack) = trapz(t(1:i), f_H_percell(1:i));
            break
        end
    end
end

% Hyield after certain experiment times
for zstack = 1:length(tmut_range)
    for i = 1:length(t)
        if t(i) >= tmut_range(zstack)
            Harea_tmut_store(ROW,COL,zstack) = trapz(t(1:i), f_H_percell(1:i));
            break
        end
    end
end

%%% Track progress

toc
X = [num2str(track/(length(prom_range)*length(z_range))*100),'%'];
disp(X);
% disp(track/(length(prom_range)*length(z_range))*100,'%')
track=track+1;

end % ROW

end % COL



%% Percentage of synthetic protein per proteome

R_frac = (R+TL_all)*n_R / mass;
C_frac = C*n_C / mass;
Q_frac = Q*n_Q / mass;
H_frac = H*n_H / mass;

R_aa = (R+TL_all)*n_R;
C_aa = C*n_C;
Q_aa = Q*n_Q;
H_aa = H*n_H;


%% Test growth rate difference via proteome

close
set(0,'defaultAxesFontSize',24);
set(gcf, 'Position',  [500, 200, 450, 320])

xx=1;

plot(t, R_frac(:,xx), 'r', 'LineWidth', 1)
hold on
plot(t, C_frac(:,xx), 'g', 'LineWidth', 1)
plot(t, Q_frac(:,xx), 'b', 'LineWidth', 1)
plot(t, H_frac(:,xx), 'k', 'LineWidth', 1)

yyaxis right
plot(R_aa(:,xx)+C_aa(:,xx)+Q_aa(:,xx)+H_aa(:,xx), 'y')
yyaxis left

ylim([0 1])
xlim([0 1000])

xxx=1;



%% Plot individual population

% close
% xx=8;
% set(gcf, 'Position',  [900, 200, 550, 500])
% 
% plot(t, subpop(:,xx), 'LineWidth', 2)
% hold on
% 
% xlim([0 80])
% xlabel('Time / h')
% ylabel('Subpop')
% % ylim([0 1e9])
% 
% axis square


%% Heat map

selection = 50;

rangemin = min(Harea_tmut_store(:,:,selection),[],'all'); % Decide on time for min
rangemax = max(Harea_tmut_store(:,:,selection),[],'all'); % For all time

% 1.5h=3 | 3h=6 | 6h=12 | 12h=24

close
set(gcf, 'Position',  [1000, 250, 300, 275])
hm = heatmap(flipud(Harea_tmut_store(:, :, selection)));

%%% Tick labels
ax = gca;
% ax.XDisplayLabels = [-12; nan(length(ax.XDisplayData)-2,1); -2];
% ax.YDisplayLabels = [5; nan(length(ax.YDisplayData)-2,1); 2.5];
ax.XDisplayLabels = nan(length(ax.XDisplayData),1);
ax.YDisplayLabels = nan(length(ax.YDisplayData),1);
S = struct(hm); 
S.XAxis.TickLabelRotation = 0;  % Don't rotate x-ticks
%%% Colour bar
caxis([rangemin rangemax]) % Range
% AX = struct(gca);
% cb = AX.Colorbar;
% cb.Ticks = [1e6 2e6 3e6];
% % cb.TickLabels = {'0','1.4e6','2.8e6'};
colorbar off
%%% General
% xlabel({'E-cell to M-cell', 'mutation probability'})
% ylabel({'E-cell growth rate','relative to M-cells'})
ax.FontSize = 24;
colormap parula

hm.GridVisible = 'off';
% hHeatmap = struct(hm).Heatmap;
% hHeatmap.GridLineStyle = ':';

% save('_Data/t_Hdrop_store_41by41.mat', 't_Hdrop_store')
% save('_Data/Harea_Hdrop_store_41by41.mat', 'Harea_Hdrop_store')
% save('_Data/Harea_tmut_store_41by41.mat', 'Harea_tmut_store')


%% Contour map for time until drop

offBlack = [38 38 38]/255;
offWhite = [217 217 217]/255;
offGrey = [100 100 100]/255;

% timeblock = [0 12 24 36 48 60 72 96];
timeblock = [0 6 12 24 48 96];

% [0.95, 0.75, 0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00001]
xx=11;

close
set(gcf, 'Position',  [1000, 250, 400, 375])
contourf(t_Hdrop_store(:, :, xx), timeblock, 'LineColor', offWhite, 'LineWidth', 1.5);
hold on
contour(t_Hdrop_store(:, :, xx),[192 192],'LineColor', offBlack, 'LineWidth', 1);

% contourf(t_Hdrop_store(:, :, xx), 'LineColor', offWhite);

xticks([]); yticks([]);
% caxis(log([3 192]));
caxis([min(timeblock), max(timeblock)]);

colormap parula


% cb = colorbar;
% cb.Ticks = log([3,12,24,96]);
% cb.TickLabels = {'3','12','24','96'};


%% Test

timeblock = [0 12 24 36 48 60 72];
rangemin = min(timeblock);
rangemax = max(timeblock);

% [0.95, 0.75, 0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00001]
xx=11;

close
set(gcf, 'Position',  [1000, 250, 300, 275])
hm = heatmap(flipud(t_Hdrop_store(:, :, xx)));

%%% Tick labels
ax = gca;
% ax.XDisplayLabels = [-12; nan(length(ax.XDisplayData)-2,1); -2];
% ax.YDisplayLabels = [5; nan(length(ax.YDisplayData)-2,1); 2.5];
ax.XDisplayLabels = nan(length(ax.XDisplayData),1);
ax.YDisplayLabels = nan(length(ax.YDisplayData),1);
S = struct(hm); 
S.XAxis.TickLabelRotation = 0;  % Don't rotate x-ticks
%%% Colour bar
caxis([rangemin rangemax]) % Range
% AX = struct(gca);
% cb = AX.Colorbar;
% cb.Ticks = [1e6 2e6 3e6];
% % cb.TickLabels = {'0','1.4e6','2.8e6'};
colorbar off
%%% General
% xlabel({'E-cell to M-cell', 'mutation probability'})
% ylabel({'E-cell growth rate','relative to M-cells'})
ax.FontSize = 24;
colormap parula

hm.GridVisible = 'off';