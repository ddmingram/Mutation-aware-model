
%%% To Do



%% Set-up

close all; clear; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',20);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
lightgrey = [0.85 0.85 0.85];
mediumgrey = [0.55 0.55 0.55];
darkgrey = [0.1 0.1 0.1];


%% Variables

%%% Global variables
N = 1e4; % Chemostat limit
num_steps = 1e8; % Steps in one simulation
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate

%%% Starting variables
t = 0; % Time
E = 1; % Engineered cell
M = 0; % Mutant cell
alpha = 0; % Burden boost
m = [3,3]; % G,H
p = [3600,3600]; % G,H

%%% Store vectors
t_vec = zeros(num_steps,1); % Time vector
EM_vec = zeros(num_steps,2); % E/M vector
p_vec = zeros(num_steps, 2);
m_vec = zeros(num_steps, 2);

%%% Record values at time 1
EM_vec(1,:) = [E,M];
p_vec(1,:) = p;
m_vec(1,:) = m;

%%% Useful other
num_rxns = 12; % No. possible reactions per cell

%%%%%
tic
%%%%%

%%
for step = 2:num_steps+1 % So we record first results as 0
    
%     step  
    
    % mRNA production
    km_G = 24/2; % molecs/h (assuming all spedies equal)
    km_H = 24/2 * (E/(E+M)); % Scale km by proportion of M!

    % mRNA decay
    dm_G = m(1) * 8/2; % molecs/h (assuming all spedies equal)
    dm_H = m(2) * 8/2;    
    
    % Protein production
    kp_G = m(1) * 24/2; % molecs/h (assuming all spedies equal)
    kp_H = m(2) * 24/2;

    % Protein decay
    dp_G = p(1) * 0.02/2; % molecs/h (assuming all spedies equal)
    dp_H = p(2) * 0.02/2;
    
    % Cell division
    if (E+M) < N
        lamE = lam0 * E*(1-z);
        lamM = lam0 * (E*z+M) * (1+alpha);
        lamEtoM = 0;
        lamMtoE = 0;
    else % Chemostat limit
        lamE = 0;
        lamM = 0;
        lamEtoM = lam0 * ((E*z+M) * E/(E+M)) * (1+alpha); % E --> M
        lamMtoE = lam0 * (E*(1-z)) * M/(E+M); % M --> E
    end
    
    rates = [(E+M)*[km_G, km_H, dm_G, dm_H, kp_G, kp_H, dp_G, dp_H],...
                    lamE, lamM, lamEtoM, lamMtoE];
    
    if E > 0 % Stop simulation once we have all M
        
        % Calculate next reaction and time from random numbers
        r1 = rand * sum(rates);
        r2 = exprnd(1/sum(rates));
        t = t + r2;
        
        % Get index of next reaction in rate_list
        % r1 < cumsum(rates): 0 when false, 1 when true
        % sum(x==0): count the zeros... +1 for first rxn
        rxn = sum(r1 < cumsum(rates)==0) + 1;     
        
        if rxn == 1
            m(1) = m(1) + 1;
        elseif rxn == 2
            m(2) = m(2) + 1;
        elseif rxn == 3
            m(1) = m(1) - 1;
        elseif rxn == 4
            m(2) = m(2) - 1;
        elseif rxn == 5
            p(1) = p(1) + 1;
        elseif rxn == 6
            p(2) = p(2) + 1;
        elseif rxn == 7
            p(1) = p(1) - 1;
        elseif rxn == 8
            p(2) = p(2) - 1;
        elseif rxn == 9
            E = E+1;
        elseif rxn == 10
            M = M+1;
        elseif rxn == 11
            E = E-1;
            M = M+1;
        elseif rxn == 12
            E = E+1;
            M = M-1;
        end        
          
    else
        ref = step;
        break
        
    end
    
    alpha = 1 - (p(end)/sum(p));

    t_vec(step) = t;
    EM_vec(step,:) = [E,M];
    m_vec(step,:) = m; 
    p_vec(step,:) = p;    
end

%%%%%
toc
%%%%%


%% Post-processing

% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec(ref:end) = [];
    EM_vec(ref:end, :) = [];
    m_vec(ref:end,:) = [];
    p_vec(ref:end,:) = [];    
end

% How alpha caries over time
alpha_vec = zeros(size(p_vec,1), 1);
for i = 1:length(alpha_vec)
    alpha_vec(i) = 1 - ( p_vec(i,end)/sum(p_vec(i,:)) );
end

% For plotting transition rate
lamEtoM = lam0 * ((EM_vec(:,1).*z + EM_vec(:,2)) .* EM_vec(:,1)./(EM_vec(:,1)+EM_vec(:,2))).*(1+alpha_vec);



%% Get evenly spaced points in time

%%% Create vector of all sequential EM transitions

EM_dynamic = zeros(num_steps,3); % t, E, M
EM_dynamic(1,:) = [0,1,0];
count = 1;
ref = EM_vec(1,:);

for i = 1:size(EM_vec, 1)    
    
    % If arrays are not equal
    if sum(EM_vec(i,:) ~= ref) > 0 % I get a TRUE/FALSE value for each array entry comparison       
        count = count + 1;
        EM_dynamic(count,:) = [t_vec(i), EM_vec(i,:)];
        ref = EM_vec(i,:); % The new reference to compare against
    end    
end

EM_dynamic(count+1:end, :) = []; % Trim


%%% Create vector of E/M at evenly spaced time points

stepsize = 0.05;
currentstep = stepsize; % Initialise

points = ceil(t_vec(end) / stepsize);
EM_tspace = zeros(points, 3); % t, E, M
EM_tspace(1,:) = [0,1,0];
count = 1; % Indexing for EM_tspace

for i = 2:size(EM_dynamic,1)
    
    time = EM_dynamic(i,1); % Next EM timepoint to check
    
    % If greater than current step, need to increment currentstep
    % for same value of E/M
    if time >= currentstep
        
        time_bound = stepsize*floor(time / stepsize); % What timepoint to we go up to?
       
        % Ensure that distance to next data point is greater than next step
        if time_bound > currentstep 
       
            for j = currentstep : stepsize: time_bound            
                 count = count + 1;            
                 EM_tspace(count,:) = [j, EM_dynamic(i-1, 2:3)];
            end

            currentstep = j + stepsize; % Add one to timebound value
        end
        
        % If next data point is before the next step, we move on to next
        % EM time point which will now reference a natural jump in EM number
       
    end   
end

%%%%%
toc
%%%%%


%% EM plot

set(gcf, 'Position',  [0, 400, 800, 640])

plot(t_vec, EM_vec(:,1), 'Color', E_yellow, 'LineWidth', 2)
hold on
plot(t_vec, EM_vec(:,2), 'Color', M_orange, 'LineWidth', 2)
legend('E', 'M')
xlabel('Time / arb. unit')
ylabel('Number of cells')
% xlim([0 8])
grid on
axis square



%% EM state space

% set(gcf, 'Position',  [0, 400, 800, 640])
% 
% plot(EM_tspace(:,2), EM_tspace(:,3), 'x', 'LineWidth', 2)
% xlabel('Number of E cells')
% ylabel('Number of M cells')
% xlim([0 N])
% grid on
% axis square



