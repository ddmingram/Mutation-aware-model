% Run the set-up file
f_setup_toggle; 

ABC_A = [95 126 196]/255;
ABC_B = [252 141 98]/255;
ABC_C = [88 178 150]/255;



%% 1. Define parts and probabilities to consider

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parts_vals_mat = [d1s1  d2s1  ...
%                   d1s2  d2s2  ...
%                     :     :      ]
%
% z_mat = [d1,->s1  d2,->s1  ...    -> = going to state x
%          d1,->s2  d2,->s2  ...  ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Record all values that I might use

%%% Set up promoter strength
part_strength = 1e4;
part_strength_int = 2.5e3;

parts_mat_all = [[part_strength; part_strength_int; 0],...      % prom_A
                 [part_strength; part_strength_int; 0],...      % RBS_A
                 [1; 0.2; 0.1],...                              % CDS_A
                 [part_strength; part_strength_int; 0],...      % prom_B
                 [part_strength; part_strength_int; 0],...      % RBS_B
                 [1; 0.2; 0.1]];                                % CDS_B

z_severe = 1e-9;
z_partial = 1e-3;

z_mat_all = [[z_severe; z_partial],...           % prom_A
             [z_severe; z_partial],...           % Etc.
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial]];

%%% Essential params
dim_vec = [1,0,0, 1,0,0];  % Which dims are active? idx=1 -> active.
d = sum(dim_vec);   % How many dims are there?
s = 3;              % How many states do we want?
n = s^d;


%% 2. Run functions to get all info needed...

%%% Get structure of the subpop network
MyPop = f_grid_structure(s,d);

%%% Order and distribute part values to params for each subpop
[z_mat, parts_mat,...
 prom_plusA_vec,...
 prom_minusA_vec,...
 RBS_plusA_vec,...
 RBS_minusA_vec,...
 CDSA_vec,...
 prom_plusB_vec,...
 prom_minusB_vec,...
 RBS_plusB_vec,...
 RBS_minusB_vec,...
 CDSB_vec] = f_PartParameters_2(dim_vec, d, s, n,...
                                parts_mat_all,...
                                z_mat_all,...
                                n_HA, n_HB);

%%% Order and distribute part values to params for each subpop
[subpops_up,...
 z_values_up,...
 subpops_dn,...
 z_values_dn] = f_popODEinfo(z_mat, MyPop, s, d, n);



%% Time variables

%%% Time testing
t_total = 1e2;

t_start = 13;
t_end = 14;
t_trials = 2;
t_switch_vec = linspace(t_start, t_end, t_trials);

%%% Store variables
InhibitorSwitch_store = zeros(1, t_trials);
t_bistable_store = zeros(1, t_trials);
% FPdistance_store = zeros(1, t_trials);
% TimeToSwitch_store = zeros(1, t_trials);

saturated_flag = false; % Check whether to terminate everything early

count_prom = 0;
% for jj = 1:t_trials
for jj = 1:1

if saturated_flag == true
    break
end

% t_switch = t_switch_vec(jj);
t_switch = 25;

tspan1 = [0 t_switch];
tspan2 = [t_switch t_total];


%% Inhibitor variables

%%% Set up inhibitor
inhibitor_min = 261;
inhibitor_max = 264;
trials_I = 4;
Inhibitor_vec = linspace(inhibitor_min,inhibitor_max,trials_I);

threshold = 50; % Percentage cells required for 'bistable'

count_ind = 0;
% for ii = 1:trials_I
for ii = 1:1

% I_set = Inhibitor_vec(ii);
I_set = 0;

%%% Collect parameters, replacing n from base_params
base_params(2) = n;
I_A = zeros(1,n); % Initialise
I_B = zeros(1,n);

params = [base_params,...
          I_A, I_B,...
          prom_plusA_vec, prom_minusA_vec,...
          RBS_plusA_vec, RBS_minusA_vec,...
          CDSA_vec,...
          prom_plusB_vec, prom_minusB_vec,...
          RBS_plusB_vec, RBS_minusB_vec,...
          CDSB_vec];

      
      
%% 4. Initialise and store variables

subpop_0 = [N, zeros(1,n-1)]; % Start with all first population (E)

e_0 = ones(1,n) * [9.86016964102649]; 
m_R_0 = ones(1,n) * [57.4983808706535];
m_C_0 = ones(1,n) * [8.01460623762431];
m_Q_0 = ones(1,n) * [972.299699264692];
m_HB_0 = ones(1,n) * 0;
TL_R_0 = ones(1,n) * [84.7428227164465];
TL_C_0 = ones(1,n) * [3.88788735587225];
TL_Q_0 = ones(1,n) * [471.662810563207];
TL_HB_0 = ones(1,n) * 0;
R_0 = ones(1,n) * [1.62351202599935];
C_0 = ones(1,n) * [1216.78308561404];
Q_0 = ones(1,n) * [147615.041146052];
HB_0 = ones(1,n) * 0;

% m_HA_0 = ones(1,n) * 1;
% TL_HA_0 = ones(1,n) * 0;
% HA_0 = ones(1,n) * 0;
m_HA_0 = ones(1,n) * [315.296521421272];
TL_HA_0 = ones(1,n) * [504.780931825327];
HA_0 = ones(1,n) * [157979.796480951];

var = [subpop_0, e_0,...
       m_R_0, m_C_0, m_Q_0, m_HA_0, m_HB_0,...
       TL_R_0, TL_C_0, TL_Q_0, TL_HA_0, TL_HB_0,...
       R_0, C_0, Q_0, HA_0, HB_0];



%% 5. Main simulations

%%% Pre signal

[t1,y1] = ode15s( @(t,y) MWCM2_odes_toggle(t, y, params, subpops_up, z_values_up, z_values_dn),...
               tspan1,...
               var,...
               Opt_1);
           
%%% Post signal

var_signal = y1(end,:); % End of tspan1 = init for tspan2
% params(length(base_params)-1) = I_set; % Change I_A - if one cell type

% Distribute inhibitor across each cell type proportionally (well-mixed)
params(length(base_params)+1 : length(base_params)+n) = I_set*y1(end,1:n)/N;

[t2,y2] = ode15s( @(t,y) MWCM2_odes_toggle(t, y, params, subpops_up, z_values_up, z_values_dn),...
               tspan2,...
               var_signal,...
               Opt_1);

%%% Stitch together
t = [t1(1:end-1); t2]; % Don't take last point to avoid double accounting
y = [y1(1:end-1,:); y2];

%%% Extract
x=1;
subpop = y(:, x:x+n-1); x=x+n;
e = y(:, x:x+n-1); x=x+n;
m_R = y(:, x:x+n-1); x=x+n;
m_C = y(:, x:x+n-1); x=x+n;
m_Q = y(:, x:x+n-1); x=x+n;
m_HA = y(:, x:x+n-1); x=x+n;
m_HB = y(:, x:x+n-1); x=x+n;
TL_R = y(:, x:x+n-1); x=x+n;
TL_C = y(:, x:x+n-1); x=x+n;
TL_Q = y(:, x:x+n-1); x=x+n;
TL_HA = y(:, x:x+n-1); x=x+n;
TL_HB = y(:, x:x+n-1); x=x+n;
R = y(:, x:x+n-1); x=x+n;
C = y(:, x:x+n-1); x=x+n;
Q = y(:, x:x+n-1); x=x+n;
HA = y(:, x:x+n-1); x=x+n;
HB = y(:, x:x+n-1);

%%% Some calcs for GR
TL_rate = (v_TL * e) ./ (K_TL + e);
gamma_HB = TL_HB .* TL_rate ./ CDSB_vec;
TL_all = TL_R + TL_C + TL_Q + TL_HA + TL_HB;
GR = TL_rate .* TL_all / mass;
d_HB = gamma_HB - (GR+p_deg).*HB;


%% Calculation - determine whether we have switched fixed points - if CONT after time x

% The following assumes that FPs always remain either side of x/y line

bistable_flag = false; % If bistable flag

HA_tot = sum(HA.*subpop,2)/N;
HB_tot = sum(HB.*subpop,2)/N;

%%% Number of cells in population where HB>HA !
% +1: Ensures there's a real difference when HA~=HB~=0
bistable_cells = sum((HB > HA+1) .* subpop, 2);

%%% Do we have >x% bistable cells over any time point?
if any( bistable_cells > N*(threshold/100))    
    bistable_flag = true;
    disp('Bistable')
    indices_above = find(bistable_cells>N*(threshold/100));
    t_bistable_store(jj) = t(indices_above(end)) - t(indices_above(1));
    InhibitorSwitch_store(jj) = Inhibitor_vec(ii);
    break
end

disp('Not bistable')
count_ind = count_ind+1;
disp(count_ind)

if I_set == Inhibitor_vec(end) % No more inhibitor is possible
    saturated_flag = true;
end





end % Inhibitor vec

count_prom = count_prom+1;
disp(count_prom)

toc

end % Time vec


%% Test plot

% close
% plot(t,bistable_cells)
% xline(t_switch)
% yline(N)
% yline(N/2,'--')
% 
% close
% plot(t, HA_tot)
% hold on
% plot(t, HB_tot)

close
set(gcf, 'Position',  [300, 50, 500, 700])
set(0,'defaultAxesFontSize',24);

subplot(2,1,1)
plot(t, HA_tot, 'linewidth', 2.5)
hold on
plot(t, HB_tot, 'color', MATLAByellow, 'linewidth', 2.5)
xline(t_switch,'--', 'linewidth', 1)
ylabel('Protein per cell')
xlim([0 100])
yyaxis right
plot(t, bistable_cells, 'k', 'linewidth', 1)
ylabel('# switched cells')
yline(N)
yline(N/2,'--', 'linewidth', 1)
ax = gca;
ax.YAxis(2).Color = 'k';
yyaxis left

subplot(2,1,2)
plot(t, subpop(:,1), 'k')
hold on
plot(t, subpop(:,2), 'color', [0.4660, 0.6740, 0.1880])
plot(t, subpop(:,5), 'color', I_orange)
plot(t, subpop(:,9), 'color', M_red)
xline(t_switch, '--', 'linewidth', 1)
xlim([0 100])
xlabel('Time')
ylabel('# cells')



%% Plot for population behaviour

close
set(gcf, 'Position',  [700, 150, 473, 271.5])
% set(gcf, 'Position',  [700, 150, 655, 376])
set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',24);

plot(t,HA_tot, 'color', ABC_A)
hold on
plot(t,HB_tot, 'color', ABC_B)
xlim([0 40])
ylim([0 2e5])
yticks([])
xticks([])
xline(t_switch, 'k--', 'LineWidth', 2)

% xlabel('Time/h')
% ylabel('Synthetic protein per cell')

yyaxis right
% a1 = area(t, subpop(:,1));
% a1.FaceColor = [1 1 1];
% a1.EdgeColor = [0 0 0];
% a1.LineWidth = 0.1;
% hold on

a2 = area(t,subpop(:,2));
a2.FaceColor = I_orange;
a2.EdgeColor = [0 0 0];
a2.LineWidth = 0.1;

% a9 = area(t,subpop(:,end));
% a9.FaceColor = M_red;
% a9.EdgeColor = [0 0 0];
% a9.LineWidth = 0.1;

alpha(0.075)

ax = gca;
ax.YAxis(2).Color = 'k';
ylim([0,1.25e9])
yticks([])
% ylabel('Number of cells in state')
yyaxis left

% exportgraphics(gca,'_Results/PaperResults/ts_s3_pop_p1e4_25e3_z1e-9_1e-3.pdf','BackgroundColor','none','ContentType','vector')



%% SI plot - protein dynamics

close
set(gcf, 'Position',  [700, 150, 500, 400])
set(0,'DefaultLineLineWidth',5);
set(0,'defaultAxesFontSize',24);

plot(t,HA_tot, 'color', ABC_A)
hold on
plot(t,HB_tot, 'color', ABC_B)
xlim([0 40])
xticks([0 40])
ylim([0 2e5])
yticks([0 2e5])
xlabel('Time/h')
ylabel('Protein per cell')
axis square



%% SI plot - cells per state

close
f = figure;
set(gcf, 'Position',  [700, 150, 800, 700])
set(0,'DefaultLineLineWidth',7);
set(0,'defaultAxesFontSize',16);

fig=1;

for i = 1:9
    
subplot(3, 3, fig)
plot(t,subpop(:,i))

xlim([0 40])
ylim([0 1e9])
xlabel('')
ylabel('')
xticks([])
yticks([])
hold on
axis square

fig=fig+1;

end

% exportgraphics(f,'_Results/PaperResults/ts_s3_allstates_p1e4_25e3_z1e-9_1e-3.pdf','BackgroundColor','none','ContentType','vector')



%% Test GR plot

% GRrel = GR./(GR(:,end));
% 
% close
% plot(t, GRrel(:,1), 'k')
% hold on
% plot(t, GRrel(:,2), 'color', [0.4660, 0.6740, 0.1880])
% plot(t, GRrel(:,5), 'color', I_orange)
% plot(t, GRrel(:,9), 'color', M_red)
% xlabel('Time/h')
% ylabel('Relative GR')
% xline(t_switch, '--', 'linewidth', 1)
% xlim([0 100])
% ylim([0 1])



%% Plot A/B conc and Phase Plane

% %%% In order to see good phase plane dots, the tspan needs to be
% %%% sufficiently short such that there are enough ODE points
% 
% close
% set(gcf, 'Position',  [800, 150, 1000, 700])
% 
% subplot(1,2,1)
% 
% plot(t, HA)
% hold on
% plot(t, HB)
% 
% xlabel('Time/h')
% ylabel('Molecule count')
% grid on
% axis square
% 
% subplot(1,2,2)
% 
% %%% Constant signal from time 0 or no signal
% % plot(HA(tspace_idx), HB(tspace_idx), 'k.')
% 
% %%% Constant signal after time x
% SignalStart = floor(tspan1(end)/multiple);
% plot(HA(tspace_idx(1:SignalStart-1)), HB(tspace_idx(1:SignalStart-1)), 'k.')
% hold on
% plot(HA(tspace_idx(SignalStart:end)), HB(tspace_idx(SignalStart:end)), 'r.')
% 
% %%% Temporary signal at time x
% % SignalIn = floor(signal(1)/multiple);
% % SignalOut = floor(signal(2)/multiple);
% % plot(HA(tspace_idx(1:SignalIn-1)), HB(tspace_idx(1:SignalIn-1)), 'k.')
% % hold on
% % plot(HA(tspace_idx(SignalIn:SignalOut)), HB(tspace_idx(SignalIn:SignalOut)), 'r.')
% % plot(HA(tspace_idx(SignalOut+1:end)), HB(tspace_idx(SignalOut+1:end)), 'k.')
% 
% xlabel('A')
% ylabel('B')
% grid on
% axis square



%% Plot the amount of inducer needed to switch

% close
% set(gcf, 'Position',  [300, 150, 1600, 600])
% 
% subplot(1,3,1)
% semilogx(prom_vec, InhibitorSwitch_store, 'kx')
% xlabel('log (Induction)')
% ylabel('Inhibitor')
% xlim([1e-1, 1e9])
% xticks([1e-1 1e1 1e3 1e5 1e7 1e9])
% xticklabels({'-1','1','3','5','7','9'})
% ylim([0 750])
% yticks([0 250 500 750])
% axis square
% % grid on
% 
% subplot(1,3,2)
% semilogx(prom_vec, TimeToSwitch_store, 'kx')
% xlabel('log (Induction)')
% ylabel('Time to switch')
% xlim([1e-1, 1e9])
% xticks([1e-1 1e1 1e3 1e5 1e7 1e9])
% xticklabels({'-1','1','3','5','7','9'})
% axis square
% 
% subplot(1,3,3)
% 
% semilogx(prom_vec, FPdistance_store, 'kx')
% xlabel('log (Induction)')
% ylabel('FP distance')
% xlim([1e-1, 1e9])
% xticks([1e-1 1e1 1e3 1e5 1e7 1e9])
% xticklabels({'-1','1','3','5','7','9'})
% ylim([0 5e7])
% axis square


%% Test steady state protein and mRNA

% close
% plot(t,m_HA)
% xlim([0 100])
% ylim([0 10])
% xlabel('Time/h')
% ylabel('mRNA')
% yyaxis right
% plot(t,HA)
% ylabel('Protein')
% ylim([0 2500])
% yticks([0 1250 2500])
% yyaxis left
% axis square

%% Legacy - PreSim steadystate

% e_0 = [9.86016020995274,9.86008187174738,9.86014336198233,23.0502435544389]; 
% m_R_0 = [57.4982913478750,57.4979306995228,57.4982190691782,12.1015042756124];
% m_C_0 = [8.01460375256224,8.01458328650765,8.01459944568171,2.50884469718596];
% m_Q_0 = [972.298981292604,972.296374152332,972.298532046440,96.5367857954715];
% m_HA_0 = [315.296505676397,0,315.296505676397,0];
% m_HB_0 = [84.7427886618350,84.7427160721661,0,0];
% TL_R_0 = [84.7427787249450,84.7427254666515,84.7427076028351,1187.57007996988];
% TL_C_0 = [3.88789140779120,3.88791039936309,3.88789240728398,69.5589263082163];
% TL_Q_0 = [471.663099697072,471.665346687761,471.663257314356,2676.52883285726];
% TL_HA_0 = [504.780787348091,0,504.780787348091,0];
% TL_HB_0 = [1.62351605113857,1.62351461657351,0,0];
% R_0 = [1.62351375890635,1.62352204229051,1.62351422754536,113.797442173062];
% C_0 = [1216.78349316197,1216.78261482363,1216.78230546059,5894.33797183080];
% Q_0 = [147615.073871017,147615.091327222,147615.069898432,226805.765550571];
% HA_0 = [157979.794645984,0,157979.794645984,0];
% HB_0 = [0,0,0,0];

% ss_store = [e(end,:);...
%          m_R(end,:);...
%          m_C(end,:);...
%          m_Q(end,:);...
%          m_HA(end,:);...
%          m_HB(end,:);...
%          TL_R(end,:);...
%          TL_C(end,:);...
%          TL_Q(end,:);...
%          TL_HA(end,:);...
%          TL_HB(end,:);...
%          R(end,:);...
%          C(end,:);...
%          Q(end,:);...
%          HA(end,:);...
%          HB(end,:)];

% % Presim for steady state values
% [t0,y0] = ode15s( @(t,y) MWCM2_odes_toggle(t, y, params, subpops_up, z_values_up, z_values_dn),...
%                [0,1e3],...
%                [subpop_0, ones(1,n*16)],...
%                Opt_1);
% x=n+1;
% e_0 = y0(end, x:x+n-1); x=x+n;
% m_R_0 = y0(end, x:x+n-1); x=x+n;
% m_C_0 = y0(end, x:x+n-1); x=x+n;
% m_Q_0 = y0(end, x:x+n-1); x=x+n;
% m_HA_0 = y0(end, x:x+n-1); x=x+n;
% m_HB_0 = zeros(1,n); x=x+n;
% TL_R_0 = y0(end, x:x+n-1); x=x+n;
% TL_C_0 = y0(end, x:x+n-1); x=x+n;
% TL_Q_0 = y0(end, x:x+n-1); x=x+n;
% TL_HA_0 = y0(end, x:x+n-1); x=x+n;
% TL_HB_0 = y0(end, x:x+n-1); x=x+n;
% R_0 = y0(end, x:x+n-1); x=x+n;
% C_0 = y0(end, x:x+n-1); x=x+n;
% Q_0 = y0(end, x:x+n-1); x=x+n;
% HA_0 = y0(end, x:x+n-1);
% HB_0 = zeros(1,n);

%% Legacy - other metrics

%%% For single cell type 

% % Get first steady state value of A
% presignal_idx = t(length(t1));
% 
% if HB(end) > HA(end)    
%     % Record inducer value
%     InhibitorSwitch_store(jj) = Inhibitor_vec(ii);    
%     % Record distance between FPs
%     FPdistance_store(jj) = sqrt((HA(length(t1))-HA(end))^2 +...
%                                (HB(length(t1))-HB(end))^2 +...
%                                (m_HA(length(t1))-m_HA(end))^2 +...
%                                (m_HB(length(t1))-m_HB(end))^2);                           
%     % Record time
%     for q = length(t1)+1:length(t)
%         % Check 1: rate of change is DECREASING (so we don't pick a value at signal start)
%         % Check 2: rate of change < 1e-7 (arbitrary)
%         if d_HB(q+1) < d_HB(q) && d_HB(q) < 1e-6
%             t_idx_FPswitch = q;
%             break
%         end
%     end
%     TimeToSwitch_store(jj) = t(t_idx_FPswitch) - t(length(t1));    
%     % Break the i-loop, move to next prom value
%     break
% end

%% Legacy - get equally spaced points in time (for phase plane)

% multiple = 0.05;
% 
% if length(t) > tspan2(end)/multiple
%     tspace_idx = ones(1, ceil(tspan2(end)/multiple));
% else % Can't have more divisions than ODE steps
%     multiple = tspan2(end)/length(t);
%     tspace_idx = ones(1, length(t));
%     disp('No. steps taken as no. ODE time steps')
% end
% 
% count = 1;
% for ii = 1:length(t)    
%     if t(ii)>multiple*count
%         while t(ii)>multiple*count % Needed if t(i) jumps by >1 'multiple'
%             count = count+1;
%             tspace_idx(count) = ii;
%         end       
%     end    
% end
