%%% Main script that simulates the mutation dynamics of a system with two
%%% synthetic genes (e.g. a toggle switch). Mutations to any of the genes'
%%% promoter regions or RBSs can be explored, and with varying degrees of
%%% mutation severity. The code has been augmented to allow for efficient
%%% testing of a toggle switch's behaviour, however these sections can be
%%% commented out. Users will most likely want to make changes in the
%%% following sections:
%%%%% 1. Construct parameters
%%%%% 2. Granularity of mutation simutations
%%%%% 4. Toggle switch inhibitor parameters
%%%%% 5./10. Exploration of toggle switch dynamics
%%%%% 7. Initial conditions for variables
%%%%% 8. ODE parameters
%%%%% 11+. Plotting
%%%%% f_params_2genes: any parameter values
%%%%% f_ODEs_2genes: mRNA repression / inhibitor behaviour
%%% For more contextual info, see "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.


%% 1. Define construct parameters

%%% Run the set-up file
f_params_2genes;

%%% Record construct-specific params. Alpha = max. TX strength (a proxy for
%%% promoter strength). Beta = mRNA-ribosome binding strength (a proxy for
%%% RBS strength). z = mutation probability. Values for 'parts' and for 'z'
%%% are stored in "part_matrix" and "z_matrix" below. _E/I/M = relating to
%%% E/I/M-states. _A/B = relating to genes A/B. By default, each gene has
%%% the same parameter values, however this can be updated within
%%% "part_matrix" and "z_matrix". This structure works for frameworks with
%%% s≤3. For frameworks with s>3, you need to add extra lines for the new
%%% I-state parameters and update "part_matrix"/"z_matrix".

alpha_E = 1e5;      % Promoter strengths in each state
alpha_I = 5e3;
alpha_M = 0;

beta_E = 60;        % RBS strengths in each state
beta_I = 6;
beta_M = 0;

z_M = 1e-6;         % Mutation probabilities for each state
z_I = 5e-2;

part_matrix = [[alpha_E; alpha_I; alpha_M],...  % prom_A
               [beta_E; beta_I; beta_M],...     % RBS_A
               [alpha_E; alpha_I; alpha_M],...  % prom_B
               [beta_E; beta_I; beta_M]];       % RBS_B

z_matrix = [[z_M; z_I],...  % prom_A
            [z_M; z_I],...  % RBS_A
            [z_M; z_I],...  % prom_B
            [z_M; z_I]];    % RBS_B        


%% 2. Define mutation granularity

%%% active_dims: vector of logicals relating to the system's parts, in
%%% order of "part_matrix" above. 1 = model mutations to that part. 0 = do
%%% not. E.g. [0,1,1,0] -> only consider mutations to the RBS of gene A and
%%% promoter of gene B.
active_dims = [1,0,1,0];
d = sum(active_dims); % The number of active dimensions

%%% s: number of states per dimension. For s>3, need to add in extra
%%% parameters, as explained in Section 1. The code currently forces each
%%% dimension to have the same number of states, but this could be changed.
s = 3;

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
 alpha_A_vec,...
 beta_A_vec,...
 alpha_B_vec,...
 beta_B_vec] = f_active_construct_params(active_dims, s, d, n,...
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


%% 4. Extra variables for toggle switch analysis

%%% Use this code to explore single simulations for the toggle switch. To
%%% remove the inhibitor's effect, set I_set to 0.

I_set = 0; % Concentration of inhibitor
t_switch = 10; % Time at which inhibitor is added
threshold = 50; % Percentage cells required for 'bistable'


%% 5. Code to explore switching behaviour of toggle switch

%%% This section is for exploring different inhibitor concentrations and
%%% times at which inhibitor is added. If you want to just explore single
%%% simulations, it can be commented-out alongside Section 10 below. To
%%% work, both this section and Section 10 must be uncommented.

%%% Variables for the strength of inhibitor
I_set = 0; % Starting value of inhibitor
I_increment = 1000; % Value to increment I_set by with each loop

%%% Variables for the time at which inhibitor is added
t_I_start = 1; % If testing t=0, need to write a non-zero number, e.g. 1e-9
t_I_end = 3;
t_I_trials = 3;
t_switch_vec = linspace(t_I_start, t_I_end, t_I_trials);

%%% Initialise a variable to track the amount of inhibitor required to
%%% cause ≥50% cells to go from more protein-A to more protein-B. Note:
%%% other metrics may be interesting to looks at, such as the distance
%%% between stable fixed points in phase space, and the time required for 
%%% 'switching' to happen following addition of inhibitor. Example code for
%%% these can be found at the end of this script.
InhibitorSwitch_store = zeros(1, t_I_trials);

count_time = 1; % Track the number of time intervals we've tested

for i = 1:t_I_trials

%%% A flag to tell the loops to move on to the next time point. Becomes
%%% 'true' when there is enough inhibitor to cause ≥50% cells to go from
%%% more protein-A to more protein-B.
flag_bistable = false; 

t_switch = t_switch_vec(i);

count_inhibitor = 1; % Track the number of inhibitor values we've tested

while flag_bistable == false % While there isn't enough inhibitor
    

%% 6. Parameters

%%% Collate all the parmeters in a list. Here, two sets of "zeros(1,n)" are
%%% added for the toggle switch simulations. These act as placeholders for
%%% the amount of two different inhibitors (I_A and I_B) in each state. In
%%% typical simulations, only one of these is used.
params = [n, base_params, alpha_A_vec, beta_A_vec, alpha_B_vec,...
          beta_B_vec, zeros(1,n), zeros(1,n)];
      
      
%% 7. Initialise variables

%%% Initial values for frameworks with [α_E,z_M,α_I,z_I] = [1e5,1e-6,5e3,
%%% 5e-2], active_dims = [1,0,1,0], and s≤3. Values are obtained by
%%% running pre-simulations with no inhibitor and allowing HA to reach a
%%% stable steady state, then storing the final values of that simulation.
%%% If you want to initialise frameworks with different starting
%%% conditions: run one simulation with I_set=0 and higher m_HA/HA than
%%% m_HB/HB -> uncomment 'ss_store_temp' below and save it to memory ->
%%% copy its values into a new 'ss_store' variable -> comment out
%%% ss_store_temp -> run a new simulation.

%%% ONLY UNCOMMENT TO INITIALISE NEW PARAMETER VALUES
% ss_store_temp = [e(end,:); m_R(end,:); m_Z(end,:); m_Q(end,:); m_HA(end,:);...
%                  m_HB(end,:); TL_R(end,:); TL_Z(end,:); TL_Q(end,:);...
%                  TL_HA(end,:); TL_HB(end,:); R(end,:); Z(end,:); Q(end,:);
%                  HA(end,:); HB(end,:)];

%%% Initial conditions for the framework types described above
if isequal(active_dims, [1,0,1,0]) % Each gene has promoter mutations
    if s == 1
        ss_store = [9.36543192517743;194.799773506685;27.5764325537007;5226.36931645169;2.32670482719244e-08;11101.6234051535;16.0138019074971;0.754158494705070;142.930410026582;6.36305732564183e-10;303.606479271620;0.0900364646377592;542.602070063601;102833.988917884;1.25363101027216e-06;218435.192898620];
    elseif s == 2
        ss_store = [9.36539089705396,9.36539089701096,9.36539089704406,22.0228139789266;194.798974447117,194.798974446255,194.798974445839,30.3764066215342;27.5764053756722,27.5764053756227,27.5764053756640,6.59647449661731;5226.36464081245,5226.36464079861,5226.36464080002,247.382421208327;11101.6124700541,-1.11155322108567e-50,11101.6124700740,-2.20640002155658e-40;2.32669839420929e-08,11101.6124700574,0,3.11266211994152e-41;16.0135121611099,16.0135121613590,16.0135121596496,1161.82945495617;0.754148087355196,0.754148087369856,0.754148087291245,71.7559044832291;142.928450753283,142.928450755939,142.928450740867,2691.00553578393;303.602289594854,-3.41791361870064e-52,303.602289569744,-2.58135856250826e-39;6.36295818819945e-10,303.602289601391,0,3.64511487557310e-40;0.0900351932663559,0.0900351932681442,0.0900351932587195,44.2837619509262;542.593496346644,542.593496355906,542.593496299941,6094.55493994412;102833.976650834,102833.976650955,102833.976651129,228559.324833847;218435.382043771,-2.09057047849404e-33,218435.382044824,1.83348442794927e-36;4.57801314708097e-07,218435.382043986,0,-3.49842694092498e-37];
    elseif s == 3
        ss_store = [9.36539089660810,9.36539089660725,9.36539089660702,9.36539089660328,11.9264673145893,11.9264673145895,9.36539089660311,11.9264673145895,22.0228139789264;194.798974522071,194.798974522002,194.798974522000,194.798974523504,129.685839976927,129.685839976928,194.798974523623,129.685839976930,30.3764066215352;27.5764053868302,27.5764053868224,27.5764053868224,27.5764053870316,21.2072374580266,21.2072374580268,27.5764053870480,21.2072374580267,6.59647449661751;5226.36464248015,5226.36464247884,5226.36464247869,5226.36464251619,1254.03014227298,1254.03014227299,5226.36464251892,1254.03014227299,247.382421208336;11101.6124745460,1.16334919597241e-09,6.31957081115203e-145,11101.6124746503,426.876760417865,-1.32440459090408e-62,11101.6124746570,426.876760427269,4.15956919774800e-89;2.32669839206620e-08,11101.6124745661,11101.6124745662,1.16334919651769e-09,9.40140121048318e-09,426.876760427269,3.31853166258782e-145,-1.13940123273215e-75,0;16.0135119550629,16.0135119551828,16.0135119551800,16.0135119514405,503.463660566740,503.463660566740,16.0135119511501,503.463660566743,1161.82945495617;0.754148077675549,0.754148077681270,0.754148077681151,0.754148077505013,26.0698073997748,26.0698073997748,0.754148077491332,26.0698073997746,71.7559044832285;142.928448906550,142.928448907639,142.928448907612,142.928448874171,1541.56449406829,1541.56449406828,142.928448871568,1541.56449406828,2691.00553578391;303.602285698055,3.18147904887521e-11,1.78210234286654e-146,303.602285630036,524.754577278597,-1.65786021285906e-62,303.602285624530,524.754577290150,4.57417266648558e-88;6.36295809803256e-10,303.602285700993,303.602285700946,3.18147904960013e-11,1.15570318543561e-08,524.754577290152,9.35817198038733e-147,-1.42627712347546e-75,0;0.0900351920730175,0.0900351920737235,0.0900351920737084,0.0900351920519863,4.35195017140591,4.35195017140589,0.0900351920502990,4.35195017140586,44.2837619509244;542.593489356426,542.593489360560,542.593489360428,542.593489233575,3347.62305187979,3347.62305187980,542.593489223681,3347.62305187978,6094.55493994407;102833.976662948,102833.976662944,102833.976662948,102833.976663006,197952.242498999,197952.242498999,102833.976663016,197952.242498999,228559.324833847;218435.382181684,2.28927972323924e-08,2.01689755972030e-12,218435.382138431,67383.7168237936,1.45929321278577e-40,218435.382135013,67383.7168252775,2.74952662389160e-67;4.57808780665907e-07,218435.382187184,218435.382186810,2.28982380942426e-08,1.48403805420080e-06,67383.7168252777,-6.78525920730846e-25,1.27857346186125e-53,0];
    else % s>3
        ss_store = ones(16,n);
        ss_store(7,:) = 0; % Set m_HB to zero to allow protein-A to increase
    end
else % Mutations to a different combination of parts
    ss_store = ones(16,n);
    ss_store(7,:) = 0; % Set m_HB to zero to allow protein-A to increase
end

%%% Assign initial values to starting variables
subpop_0 = [N, zeros(1,n-1)]; % Assume population always starts as all-E-state
e_0 = ones(1,n).*ss_store(1,:); 
m_R_0 = ones(1,n).*ss_store(2,:);
m_Z_0 = ones(1,n).*ss_store(3,:);
m_Q_0 = ones(1,n).*ss_store(4,:);
m_HA_0 = ones(1,n).*ss_store(5,:);
m_HB_0 = ones(1,n).*ss_store(6,:);
TL_R_0 = ones(1,n).*ss_store(7,:);
TL_Z_0 = ones(1,n).*ss_store(8,:);
TL_Q_0 = ones(1,n).*ss_store(9,:);
TL_HA_0 = ones(1,n).*ss_store(10,:);
TL_HB_0 = ones(1,n).*ss_store(11,:);
R_0 = ones(1,n).*ss_store(12,:);
Z_0 = ones(1,n).*ss_store(13,:);
Q_0 = ones(1,n).*ss_store(14,:);
HA_0 = ones(1,n).*ss_store(15,:);
HB_0 = ones(1,n).*ss_store(16,:);

%%% Save initial values as a variable
var1 = [subpop_0, e_0, m_R_0, m_Z_0, m_Q_0, m_HA_0, m_HB_0, TL_R_0,...
       TL_Z_0, TL_Q_0, TL_HA_0, TL_HB_0, R_0, Z_0, Q_0, HA_0, HB_0];


%% 8. Call ODE solver

%%% Call the 'stiff' ODE solver, ode15s. For toggle switch simulations, we
%%% simulate the effects of inhibitor appearing instantly at a t=t_switch.
%%% For this to occur, we simulate two sets of ODEs: one for [0, tswitch],
%%% and another for [t_switch, t_end], and then join them together. The
%%% initial conditions of the second set are the final variable quantities
%%% of the first set. Note: to change abs/rel toleraces, use
%%% "odeset('RelTol',x,'AbsTol',y)" as a final option.

%%% Simulations up until inhibitor is added
tspan1 = [0, t_switch];
[t1,y1] = ode15s(@(t,y) f_ODEs_2genes(t, y, params, states_up, z_values_up, z_values_dn),...
               tspan1,...
               var1);

%%% Initial conditions of second set = final values of first set
var_2 = y1(end,:);

%%% Distribute inhibitor across each cell type proportionally
I_per_state = I_set * y1(end,1:n) / N;
params(length(params)-2*n+1 : length(params)-n) = I_per_state;

%%% Simulations after inhibitor is added
t_end = 100;
tspan2 = [t_switch, t_end];
[t2,y2] = ode15s( @(t,y) f_ODEs_2genes(t, y, params, states_up, z_values_up, z_values_dn),...
               tspan2,...
               var_2);

%%% Combine simulations. '-1' is used to that the variables at the point
%%% of inhibitor entry aren't counted twice
t = [t1(1:end-1); t2];
y = [y1(1:end-1,:); y2];


%% 9. Extract variables and calculate any terms required

x = 1; % Initialise for automatic numbering
subpop = y(:, x:x+n-1); x=x+n;
e = y(:, x:x+n-1); x=x+n;
m_R = y(:, x:x+n-1); x=x+n;
m_Z = y(:, x:x+n-1); x=x+n;
m_Q = y(:, x:x+n-1); x=x+n;
m_HA = y(:, x:x+n-1); x=x+n;
m_HB = y(:, x:x+n-1); x=x+n;
TL_R = y(:, x:x+n-1); x=x+n;
TL_Z = y(:, x:x+n-1); x=x+n;
TL_Q = y(:, x:x+n-1); x=x+n;
TL_HA = y(:, x:x+n-1); x=x+n;
TL_HB = y(:, x:x+n-1); x=x+n;
R = y(:, x:x+n-1); x=x+n;
Z = y(:, x:x+n-1); x=x+n;
Q = y(:, x:x+n-1); x=x+n;
HA = y(:, x:x+n-1); x=x+n;
HB = y(:, x:x+n-1);

%%% Calculate growth rate (note: to analyse more rates, copy the relevant
%%% code from the ODE script and insert below)
TL_rate = (v_TL * e) ./ (K_TL + e);
TL_all = TL_R + TL_Z + TL_Q + TL_HA + TL_HB;
GR = TL_rate .* TL_all / mass;

%%% Per-cell average of HA and HB in the popultion
HA_avg_percell = sum(HA.*subpop,2)/N;
HB_avg_percell = sum(HB.*subpop,2)/N;

%%% Number of cells in population where HB>HA. '+1' ensures there's a real
%%% difference when HA~=HB~=0
bistable_cells = sum((HB > HA+1) .* subpop, 2);


%% 10. Toggle switch metric calculations

%%% This section is linked to the set-up in Section 5 above, and determines
%%% whether a certain number of cells have switched from more-A to more-B.
%%% To test this metric, this block should be uncommented (until the 'end'
%%% for the inhibitor loop) alongside Section 5.

%%% Do we have >x% bistable cells over any time point?
if any(bistable_cells > N*(threshold/100))
    disp('Bistable')
    flag_bistable = true; % To exit the while loop
    InhibitorSwitch_store(i) = I_set;
else
    disp('Not bistable')
    disp(count_inhibitor)
    count_inhibitor = count_inhibitor+1;
    I_set = I_set + I_increment; % Test higher values
end


end % While loop that increments inhibitor concentration until 'bistable'

disp(count_time)
count_time = count_time+1;

end % Loop to add inhibitor at different time intervals


%% 11+. Plotting

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

%%% Plot styles
set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20);


%% Plot: Number of cells where HA>HB over time

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% 
% plot(t, bistable_cells)
% xline(t_switch)
% xlim([0 32])
% xlabel('Time/h')
% yline(N)
% yline(N/2,'--')
% ylim([0 N])
% ylabel('No. switched cells')


%% Plot: No. cells in specified state

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% xx=1; % Choose state number in framework. E.g. 1 always = E-state
% 
% plot(t, subpop(:,xx))
% xlim([0 32])
% xlabel('Time/h')
% ylim([0 1e9])
% ylabel('No. cells')


%% Plot: Protein conc. per cell in specified state

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% xx=1; % Choose state number in framework. E.g. 1 always = E-state
% 
% plot(t, HA(:,xx), 'Color', ABC_A)
% hold on
% plot(t, HB(:,xx), 'Color', ABC_B)
% xlim([0 32])
% xlabel('Time/h')
% ylim([0 3e5])
% ylabel('Protein')


%% Plot: Avg protein conc. per cell (left), no. cells in specified state (right)

% close
% set(gcf, 'Position',  [500, 200, 450, 325])
% 
% plot(t, HA_avg_percell, 'color', ABC_A)
% hold on
% plot(t, HB_avg_percell, 'color', ABC_B)
% xlim([0 32])
% xlabel('Time/h')
% ylim([0 3e5])
% ylabel('Avg protein')
% 
% yyaxis right
% area_E = area(t, subpop(:,1));
% area_E.FaceColor = [1 1 1];
% area_E.EdgeColor = [0 0 0];
% area_E.LineWidth = 0.1;
% 
% %%% Uncomment to display 'intermediate' area for when s>2
% % area_I = area(t, subpop(:,5));
% % area_I.FaceColor = [224 135 0]/255;
% % area_I.EdgeColor = [0 0 0];
% % area_I.LineWidth = 0.1;
% 
% area_M = area(t,subpop(:,end));
% area_M.FaceColor = M_red;
% area_M.EdgeColor = [0 0 0];
% area_M.LineWidth = 0.1;
% 
% alpha(0.05)
% 
% ax = gca;
% ax.YAxis(2).Color = 'k';
% ylim([0,1.15e9])
% ylabel('No. cells')
% yyaxis left


%% Legacy - other metrics

%%% This is legacy code and will need adjusting, but is included as an
%%% initial guide.

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