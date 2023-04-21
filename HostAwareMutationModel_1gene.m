%%% Main script that simulates the mutation dynamics of a system with one
%%% synthetic gene. Mutations to the gene's promoter region or RBS can be explored, and
%%% with varying degrees of mutation severity. The code has been augmented
%%% to allow analysis of protein yield and cell viability when alpha_E and
%%% z are varied. Users will most likely want to make changes in the
%%% following sections:
%%%%% 1. Construct parameters
%%%%% 2. Granularity of mutation simutations
%%%%% 4./8. Code to loop simulations over z/alpha
%%%%% 6. Initial conditions for variables, ODE parameters
%%%%% 9+. Plotting
%%%%% f_params_1gene: any parameter values
%%%%% f_ODEs_1gene: mRNA repression / inhibitor behaviour
%%% For more contextual info, see "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.


%% 1. Define construct parameters

%%% Run the set-up file
f_params_1gene;

%%% Record construct-specific params. Alpha = max. TX strength (a proxy for
%%% promoter strength). Beta = mRNA-ribosome binding strength (a proxy for
%%% RBS strength). z = mutation probability. Values for 'parts' and for 'z'
%%% are stored in "part_matrix" and "z_matrix" below. _E/I/M = relating to
%%% E/I/M-states. This structure works for frameworks with s≤3. For
%%% frameworks with s>3, you need to add extra lines for the new I-state
%%% parameters and update "part_matrix"/"z_matrix".

alpha_E = 1e4;      % Promoter strengths in each state
alpha_I = 1e3;
alpha_M = 0;

beta_E = 60;        % RBS strengths in each state
beta_I = 6;
beta_M = 0;

z_M = 1e-6;         % Mutation probabilities for each state
z_I = 1e-2;

part_matrix = [[alpha_E; alpha_I; alpha_M],...  % prom
               [beta_E; beta_I; beta_M]];       % RBS

z_matrix = [[z_M; z_I],...  % prom
            [z_M; z_I]];    % RBS


%% 2. Define mutation granularity

%%% active_dims: vector of logicals relating to the system's parts, in
%%% order of "part_matrix" above. 1 = model mutations to that part. 0 = do
%%% not. E.g. [1,0] -> only consider mutations to the promoter.
active_dims = [1,0];
d = sum(active_dims); % The number of active dimensions

%%% s: number of states per dimension. For s>3, need to add in extra
%%% parameters, as explained in Section 1. The code currently forces each
%%% dimension to have the same number of states, but this could be changed.
s = 2;

n = s^d; % Number of mutation states in the framework


%% 3. Create the framework structure and extract relevant connections

%%% f_grid_structure: function that labels each state systematically with a
%%% coordinate ('coord') of length 's', then creates the cell array
%%% state_connections that collects info on how each state is connected to
%%% one another. This is required for automatically assigning 'part' and
%%% 'z' values to each state in the correct order. state_connections is
%%% structured as follows: Rows: states in order. Columns: (1) index of
%%% state. (2) coord of state. (3) List of coords of upstream states. (4)
%%% List of index pairs that each specify a value from the z_matrix for th
%%% mutation probability of an upstream state. (5-6) Same as (3-4) but for
%%% downstream states.
state_connections = f_allocate_coords(s, d, n);

%%% Obtain the 'active' parameters for part strengths and z_values, and
%%% distribute these across the states in the correct order
[active_z_matrix,...
 alpha_vec,...
 beta_vec] = f_active_construct_params(active_dims, s, d, n,...
                                         part_matrix,...
                                         z_matrix);

%%% The ODE for each state involves terms regarding (i) number of cells in
%%% connected upstream states, (ii) mutation probabilities from connected
%%% upstream states, (iii) mutation probabilities to downstream states. The
%%% cell array 'state_connections' has this info in an over-extensive and
%%% coded format. Hence, f_extract_state_connections transfers therequired
%%% info into defined variables and converts 'coordinates' to numbered
%%% indices and 'z_matrix index pairs' into probabilities.
[states_up,...
 z_values_up,...
 z_values_dn] = f_extract_state_connections(state_connections, active_z_matrix, s, d, n);


%% 4. Initialise metrics to explore protein production

%%% Uncomment this section and Section 8 below if you want to explore the
%%% effects of varying z/alpha for s2 frameworks. Can also explore effects
%%% for s3 frameworks, but code adjustments are required (see below).

% %%% Initialise the range of values for the two key variables: z ('x-axis',
% %%% columns of heat map) and alpha ('y-axis', rows of heat map). Increase
% %%% 'divisions' to have more granularity.
% divisions = 5; % Suggestion: 11 for heatmaps, 21 for contour maps
% z_range = logspace(-12, -3, divisions);
% alpha_range = logspace(2.5, 5, divisions);
% 
% %%% For each z/alpha condition, record (i) the yield of H when stopping the
% %%% experiment at different time points (yield_t_stop), and the time taken
% %%% for H to fall to certain %'s of its max (t_maxH_drop).
% t_stop_range = [0.1, linspace(1,48,48), 200]; % Final value must be within tspan
% yield_t_stop_store = zeros(length(z_range), length(alpha_range), length(t_stop_range)); % Yield store
% maxH_drop_range = [0.95, 0.90, 0.75, 0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00001];
% t_maxH_drop_store = zeros(length(z_range), length(alpha_range), length(maxH_drop_range)); % Time store
% 
% count = 1; % To track progress for large iterations
% 
% for COL = 1:length(z_range)
% 
% %%% Adjust z-values based on the loops. This works for s=2, but for s=3
% %%% need to re-call 'f_extract_state_connections' function here.
% z_values_up{2} = z_range(COL); % Need to adjust z based on loop
% z_values_dn{1} = z_range(COL); % Need to adjust z based on loop
% 
% for ROW = 1:length(alpha_range)
% 
% %%% Adjust alpha based on loop
% alpha_vec(1) = alpha_range(ROW);

    
%% 5. Collate all the parmeters in a list

params = [n, base_params, alpha_vec, beta_vec];


%% 6 . Initialise variables and call solver

%%% Initial values for frameworks with active_dims = [1,0] and s≤3. Values
%%% are chosen such that simulations start with no synthetic protein
%%% expression. Values are obtained by running pre-simulations with alpha
%%% and beta parameters set to zero, then storing the final values of that
%%% simulation. If you want to initialise frameworks with different
%%% starting conditions: run one simulation -> uncomment 'ss_store_temp'
%%% below and save it to memory -> copy its values into a new 'ss_store'
%%% variable -> comment out ss_store_temp -> run a new simulation.

%%% ONLY UNCOMMENT TO INITIALISE NEW PARAMETER VALUES
% ss_store_temp = [e(end,:); m_R(end,:); m_Z(end,:); m_Q(end,:); m_H(end,:);...
%                  TL_R(end,:); TL_Z(end,:); TL_Q(end,:); TL_H(end,:);...
%                  R(end,:); Z(end,:); Q(end,:); H(end,:)];

%%% Initial conditions for the framework types described above
if isequal(active_dims, [1,0]) % Promoter mutations
    if s == 1
        ss_store = [22.0228139789267;30.3764066215362;6.59647449661767;247.382421208341;2.64492288189647e-34;1161.82945495617;71.7559044832282;2691.00553578390;-2.29849403825534e-35;44.2837619509233;6094.55493994408;228559.324833847;-1.81302615717280e-21];
    elseif s == 2
        ss_store = [22.0228139789269,22.0228139789265;30.3764066215354,30.3764066215341;6.59647449661753,6.59647449661728;247.382421208335,247.382421208326;0,0;1161.82945495618,1161.82945495617;71.7559044832286,71.7559044832290;2691.00553578391,2691.00553578392;0,0;44.2837619509246,44.2837619509263;6094.55493994411,6094.55493994411;228559.324833847,228559.324833847;0,0];
    elseif s == 3
        ss_store = [22.0228139789230,22.0228139789265,22.0228139787810;30.3764066216332,30.3764066215356,30.3764066281090;6.59647449663543,6.59647449661757,6.59647449781615;247.382421209070,247.382421208337,247.382421257365;0,0,0;1161.82945495602,1161.82945495617,1161.82945494858;71.7559044831840,71.7559044832284,71.7559044803003;2691.00553578293,2691.00553578391,2691.00553571833;0,0,0;44.2837619507754,44.2837619509240,44.2837619410083;6094.55493994203,6094.55493994409,6094.55493981599;228559.324833830,228559.324833847,228559.324832706;0,0,0];
    else % s>3
        ss_store = ones(13,n);
    end
else % Mutations to a different combination of parts
    ss_store = ones(13,n);
end

%%% Assign initial values to starting variables
subpop_0 = [N, zeros(1,n-1)]; % Assume population always starts as all-E-state
e_0 = ones(1,n).*ss_store(1,:); 
m_R_0 = ones(1,n).*ss_store(2,:);
m_Z_0 = ones(1,n).*ss_store(3,:);
m_Q_0 = ones(1,n).*ss_store(4,:);
m_H_0 = ones(1,n).*ss_store(5,:);
TL_R_0 = ones(1,n).*ss_store(6,:);
TL_Z_0 = ones(1,n).*ss_store(7,:);
TL_Q_0 = ones(1,n).*ss_store(8,:);
TL_H_0 = ones(1,n).*ss_store(9,:);
R_0 = ones(1,n).*ss_store(10,:);
Z_0 = ones(1,n).*ss_store(11,:);
Q_0 = ones(1,n).*ss_store(12,:);
H_0 = ones(1,n).*ss_store(13,:);

%%% Save initial values as a variable
var = [subpop_0, e_0, m_R_0, m_Z_0, m_Q_0, m_H_0, TL_R_0, TL_Z_0,...
       TL_Q_0, TL_H_0, R_0, Z_0, Q_0, H_0];

%%% Call the 'stiff' ODE solver, ode15s. Note: to change abs/rel
%%% toleraces, use "odeset('RelTol',x,'AbsTol',y)" as a final option
tspan = [0, 200]; % Simulation time range
[t,y] = ode15s(@(t,y) f_ODEs_1gene(t, y, params, states_up, z_values_up, z_values_dn),...
               tspan,...
               var);

           
%% 7. Extract variables and calculate metrics

x = 1; % Initialise for automatic numbering
subpop = y(:, x:x+n-1); x=x+n;
e = y(:, x:x+n-1); x=x+n;
m_R = y(:, x:x+n-1); x=x+n;
m_Z = y(:, x:x+n-1); x=x+n;
m_Q = y(:, x:x+n-1); x=x+n;
m_H = y(:, x:x+n-1); x=x+n;
TL_R = y(:, x:x+n-1); x=x+n;
TL_Z = y(:, x:x+n-1); x=x+n;
TL_Q = y(:, x:x+n-1); x=x+n;
TL_H = y(:, x:x+n-1); x=x+n;
R = y(:, x:x+n-1); x=x+n;
Z = y(:, x:x+n-1); x=x+n;
Q = y(:, x:x+n-1); x=x+n;
H = y(:, x:x+n-1);

%%% Calculate growth rate (note: to analyse more rates, copy the relevant
%%% code from the ODE script and insert below)
TL_rate = (v_TL * e) ./ (K_TL + e);
TL_all = TL_R + TL_Z + TL_Q + TL_H;
GR = TL_rate .* TL_all / mass;

%%% Calculate metrics to do with synthetic protein production. The 'rate' 
%%% of H-production is used to calculate the 'yield' over different time
%%% spans (trapz function below).
Hrate_avg_percell = sum(GR .* H .* subpop, 2) ./ (N * log(2));
[Hmax, Hmax_t_idx] = max(Hrate_avg_percell);


%% 8. Exploration of protein production metrics

%%% This section goes in-line with the set-up from Section 4. They must be
%%% commented/uncommented together.

% %%% H-protein yield after different times
% for zstack = 1:length(t_stop_range)
%     for i = 1:length(t)
%         if t(i) >= t_stop_range(zstack)
%             yield_t_stop_store(ROW,COL,zstack) = trapz(t(1:i), Hrate_avg_percell(1:i));
%             break % Go to the next 'time'
%         end
%     end
% end
% 
% %%% Time for H-protein to drop to % of max
% for zstack = 1:length(maxH_drop_range)
%     for i = 1:length(Hrate_avg_percell)
%         if (i > Hmax_t_idx) && (Hrate_avg_percell(i) <= maxH_drop_range(zstack)*Hmax)
%             t_maxH_drop_store(ROW,COL,zstack) = t(i);
%             break
%         end
%     end
% end
% 
% disp([num2str(count/(length(alpha_range)*length(z_range))*100),'%']); % Display progress
% count=count+1;
% 
% end % ROW (alpha_E)
% 
% end % COL (z)


%% 9+. Plotting

%%% The first section here contains a list of definitions and commands that
%%% are useful for all general plots. These can be modified or added to as
%%% appropriate. The sections after this each define a different plot that
%%% can be useful for analysis. Uncomment any section to run.

%%% Colours
ABC_A = [95 126 196]/255;
ABC_B = [252 141 98]/255;
ABC_C = [88 178 150]/255;
GRgreen = [26, 173, 27]/255;
E_yellow = [224 199 5]/255;
I_orange = [224 135 0]/255;
M_red = [255, 51, 51]/255;
offBlack = [38 38 38]/255; % For contour map
offWhite = [217 217 217]/255; % For contour map

%%% Plot styles
set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20);


%% Plot: No. cells in specified state

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% xx=1; % Choose state number in framework. E.g. 1 always = E-state
% 
% plot(t, subpop(:,xx))
% xlim([0 40])
% xlabel('Time/h')
% ylim([0 1e9])
% ylabel('No. cells')


%% Plot: Protein conc. per cell in specified state

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% xx=1; % Choose state number in framework. E.g. 1 always = E-state
% 
% plot(t, H(:,xx), 'Color', ABC_A)
% xlim([0 40])
% xlabel('Time/h')
% ylim([0 3e5])
% ylabel('Protein')


%% Heat map for H-protein yield with varying z/alpha

%%% Heat map. Required simulating loop conditions above.

% close
% set(gcf, 'Position',  [1000, 250, 500, 450])
% 
% %%% Choose which t-value to plot (see t_stop_range above).
% %%% 1='near 0', 2=1h, 3=2h etc.
% t_choice = 10;
% 
% %%% Fix the range of each heatmap to be the min/max value of that
% %%% particular simulation. This is so each heatmap shows the full gradation
% %%% of colours.
% rangemin = min(yield_t_stop_store(:,:,t_choice),[],'all'); % Decide on time for min
% rangemax = max(yield_t_stop_store(:,:,t_choice),[],'all'); % For all time
% 
% hm = heatmap(flipud(yield_t_stop_store(:, :, t_choice))); % Call heatmap
% ax = gca;
% ax.XDisplayLabels = [-12; nan(length(ax.XDisplayData)-2,1); -3]; % Normal axis labels
% ax.YDisplayLabels = [5; nan(length(ax.YDisplayData)-2,1); 2.5];
% % ax.XDisplayLabels = nan(length(ax.XDisplayData),1); % No axis labels
% % ax.YDisplayLabels = nan(length(ax.YDisplayData),1);
% caxis([rangemin rangemax]) % Range for colourbar
% xlabel('10^{alpha_E}')
% ylabel('10^{z_M}')
% ax.FontSize = 16;
% colormap parula
% hm.GridVisible = 'off';


%% Contour map for 'time until H drops to % of max' with varying z/alpha

%%% Contour map. Required simulating loop conditions above.

% close
% set(gcf, 'Position',  [1000, 250, 400, 375])
% 
% %%% Choose time-contours to plot - '192' is a separate addition (see below)
% t_contours = [0 6 12 24 48 96]; 
% 
% %%% Choose which t-value to plot (see maxH_drop_range above).
% %%% 1=95%, 2=90%, 3=75% etc.
% H_choice = 2;
% 
% %%% Base contour map
% contourf(t_maxH_drop_store(:, :, H_choice), t_contours, 'LineColor', offWhite, 'LineWidth', 1.5);
% 
% %%% Add a single contour - useful for specifying a single different colour,
% %%% e.g. when we want bigger contrast between two lighter regions.
% hold on
% contour(t_maxH_drop_store(:, :, H_choice),[192 192],'LineColor', offBlack, 'LineWidth', 1);
% 
% xticks([]); yticks([]);
% caxis([min(t_contours), max(t_contours)]);
% colormap parula