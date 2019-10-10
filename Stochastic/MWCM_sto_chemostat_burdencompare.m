%% Set-up

close all; clear; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',20);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
lightgrey = [0.85 0.85 0.85];
mediumgrey = [0.55 0.55 0.55];
darkgrey = [0.1 0.1 0.1];

%%% Global variables
N = 1e4; % Chemostat limit
num_steps = 5e8; % Steps in one simulation
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate
num_rxns = 10; % No. possible reactions per cell


%% Sim W/O burden

%%% Starting variables
t = 0; % Time
E = N; % Engineered cell
M = N-E; % Mutant cell
m = [3,3]; % G,H
p = [3600,3600]; % G,H

%%% Store vectors
t_vec_x = zeros(num_steps,1); % Time vector
EM_vec_x = zeros(num_steps, 2); % E/M vector
p_vec_x = zeros(num_steps, 2);
m_vec_x = zeros(num_steps, 2);

%%% Record values at time 1
EM_vec_x(1,:) = [E,M];
p_vec_x(1,:) = p;
m_vec_x(1,:) = m;

%%%%%
tic
%%%%%

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
    lamEtoM = lam0 * ((E*z+M) * E/(E+M)); % E --> M
    lamMtoE = lam0 * (E*(1-z)) * M/(E+M); % M --> E
    
    rates = [(E+M)*[km_G, km_H, dm_G, dm_H, kp_G, kp_H, dp_G, dp_H],...
                    lamEtoM, lamMtoE];
    
    if E > 0 % Stop simulation once we have all E
        
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
            E = E-1;
            M = M+1;
        elseif rxn == 10
            E = E+1;
            M = M-1;          
        end        
          
    else
        ref = step;
        break
        
    end
    
    t_vec_x(step) = t;
    EM_vec_x(step,:) = [E,M];
    m_vec_x(step,:) = m; 
    p_vec_x(step,:) = p;          
end


% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec_x(ref:end) = [];
    EM_vec_x(ref:end, :) = [];
    m_vec_x(ref:end,:) = [];
    p_vec_x(ref:end,:) = [];   
end

alpha_vec_x = zeros(size(p_vec_x,1), 1);
for i = 1:length(alpha_vec_x)
    alpha_vec_x(i) = 1 - ( p_vec_x(i,end)/sum(p_vec_x(i,:)) );
end

%%%%%
toc
%%%%%


%% Sim with burden

%%% Starting variables
t = 0;
E = N; % Engineered cell
M = N-E; % Mutant cell
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
alpha = 1 - p_vec(1,end) / sum(p_vec(1,:)); % Burden boost


for step = 2:num_steps+1 % So we record first results as 0
    
%     step  
    
    % Cell division
    lamEtoM = lam0 * ((E*z+M) * E/(E+M)) * (1+alpha); % E --> M
    lamMtoE = lam0 * (E*(1-z)) * M/(E+M); % M --> E
    
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
    
    rates = [(E+M)*[km_G, km_H, dm_G, dm_H, kp_G, kp_H, dp_G, dp_H],...
                    lamEtoM, lamMtoE];
    
    if E > 0 % Stop simulation once we have all E
        
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
            E = E-1;
            M = M+1;
        elseif rxn == 10
            E = E+1;
            M = M-1;          
        end        
          
    else
        ref = step;
        break
        
    end
    
    alpha = 1 - p(end)/sum(p);
    
    t_vec(step) = t;
    EM_vec(step,:) = [E,M];
    m_vec(step,:) = m; 
    p_vec(step,:) = p;        
end

% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec(ref:end) = [];
    EM_vec(ref:end, :) = [];
    m_vec(ref:end,:) = [];
    p_vec(ref:end,:) = [];
end

alpha_vec = zeros(size(p_vec,1), 1);
for i = 1:length(alpha_vec)
    alpha_vec(i) = 1 - ( p_vec(i,end)/sum(p_vec(i,:)) );
end

lamEtoM = lam0 * ((EM_vec(:,1).*z + EM_vec(:,2)) .* EM_vec(:,1)./(EM_vec(:,1)+EM_vec(:,2))).*(1+alpha_vec);

%%%%%
toc
%%%%%


%% Plot without burden

close all
set(gcf, 'Position',  [0, 400, 1400, 640])

subplot(1,2,1)

plot(t_vec_x, EM_vec_x(:,1), '--', 'Color', E_yellow)
hold on
plot(t_vec_x, EM_vec_x(:,2), '--', 'Color', M_orange)
legend('E (no burden)', 'M (no burden)', 'Location', 'East')
xlabel('Time / h')
ylabel('Number of cells')
grid on

subplot(1,2,2)

plot(t_vec_x, 1-alpha_vec_x, 'Color', lightgrey, 'LineWidth', 0.5)
legend('No burden', 'Location', 'East')
xlabel('Time / h')
ylabel('H fraction')
grid on
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))



%% Plot with burden

close all
set(gcf, 'Position',  [0, 400, 1400, 640])

subplot(1,2,1)

plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
hold on
plot(t_vec, EM_vec(:,2), 'Color', M_orange)
legend('E (no burden)', 'M (no burden)', 'Location', 'East')
xlabel('Time / h')
ylabel('Number of cells')
grid on

subplot(1,2,2)

plot(t_vec, 1-alpha_vec, 'Color', mediumgrey, 'LineWidth', 0.5)
legend('No burden', 'Location', 'East')
xlabel('Time / h')
ylabel('H fraction')
grid on
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))



%% Plot BOTH

close all
set(gcf, 'Position',  [0, 400, 1400, 640])

subplot(1,2,1)

plot(t_vec_x, EM_vec_x(:,1), '--', 'Color', E_yellow)
hold on
plot(t_vec_x, EM_vec_x(:,2), '--', 'Color', M_orange)
plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
plot(t_vec, EM_vec(:,2), 'Color', M_orange)
legend('E (no burden)', 'M (no burden)','E (burden)','M (burden)', 'Location', 'East')
reorderLegend([2,1,4,3])
% xlim([0 8])
xlabel('Time / h')
ylabel('Number of cells')
grid on

subplot(1,2,2)

plot(t_vec_x, 1-alpha_vec_x, 'Color', lightgrey, 'LineWidth', 0.5)
hold on
plot(t_vec, 1-alpha_vec, 'Color', mediumgrey, 'LineWidth', 0.5)
legend('No burden', 'Burden', 'Location', 'East')
reorderLegend([1 2])
% xlim([0 8])
xlabel('Time / h')
ylabel('H fraction')
grid on
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))