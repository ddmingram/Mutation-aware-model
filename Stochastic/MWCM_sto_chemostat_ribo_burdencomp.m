
%%% To Do


%% Set-up

close all; clear; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',25);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
lightgrey = [0.85 0.85 0.85];
mediumgrey = [0.55 0.55 0.55];
darkgrey = [0.1 0.1 0.1];
R_blue = [0 114 189]/255;
G_green = [119 172 48]/255;
H_yellow = [237 177 32]/255;


%% Global Variables

%%% Global variables
N = 1e4; % Chemostat limit
num_steps = 1e8; % Steps in one simulation
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate

%%% Useful other
num_rxns = 14; % No. possible reactions per cell

%%% NEW TESTING
km_tot = 24;
dm_tot = 8;
kp_tot = 24;
dp_tot = 0.02;
rates_vec = zeros(num_steps,14);



%% Sim LESS burden 

%%% Starting variables
t = 0; % Time
E = N; % Engineered cell
M = E-N; % Mutant cell
m = [3,3,3]; % G,H
p = [1200,1200,1200]; % R,G,H

%%% Store vectors
t_vec_x = zeros(num_steps, 1); % Time vector
m_vec_x = zeros(num_steps, 3);
p_vec_x = zeros(num_steps, 3);
EM_vec_x = zeros(num_steps, 2); % E/M vector

%%% Record values at time 1
EM_vec_x(1,:) = [E,M];
p_vec_x(1,:) = p;
m_vec_x(1,:) = m;
alpha_x = 0.1 * (1 - p_vec_x(1,end) / sum(p_vec_x(1,:))); % Lower burden

%%%%%
tic
%%%%%

%%%
for step = 2:num_steps+1 % So we record first results as 0
       
%     step   

    % Assume all species equal / rates are molecs/h for a single
    % transcript - dynamcis/trends may vary if we scale to a whole cell
    
    %%% mRNA production
    gamma_km = (51-25*(E/(E+M))) / 26;
    km_R = km_tot/51 * gamma_km;
    km_G = km_tot*25/51 * gamma_km;
    km_H = km_tot*25/51 * (E/(E+M)); % Scale km by proportion of M!
    
    %%% mRNA decay   
    dm_R = m(1) * dm_tot/51;
    dm_G = m(2) * dm_tot*25/51;
    dm_H = m(3) * dm_tot*25/51;    
    
    %%% Protein production
    if p(1) > 0 && m(1)+m(2)+m(3) > 0
        gamma_kp = 51 / (p(1)*(m(1) + 25*m(2) + 25*m(3)));        
    else
        gamma_kp = 0;
    end
    kp_R = m(1) * p(1) * kp_tot/51 * gamma_kp;
    kp_G = m(2) * p(1) * kp_tot*25/51 * gamma_kp;
    kp_H = m(3) * p(1) * kp_tot*25/51 * gamma_kp;

    %%% Protein decay 
    dp_R = p(1) * dp_tot/51;
    dp_G = p(2) * dp_tot*25/51;
    dp_H = p(3) * dp_tot*25/51;
    
    %%% Cell division
    lamEtoM = lam0 * (((E*z/(1+alpha_x)) + M) * E/(E+M)) ;%* (1+alpha_x); % E --> M
    lamMtoE = lam0 * (E*(1-z)) * M/(E+M) / (1+alpha_x); % M --> E
    
    % Record
    rates = [(E+M)*[km_R, km_G, km_H, dm_R, dm_G, dm_H, ...
                    kp_R, kp_G, kp_H, dp_R, dp_G, dp_H],...
                    lamEtoM, lamMtoE];                 
    rates_vec(step-1, :) = [rates(1:num_rxns-2)/(E+M), rates(num_rxns-1:num_rxns)];
    
    if E > 0 
        
        % Calculate next reaction and time from random numbers
        r1 = rand * sum(rates);
        r2 = -(1/sum(rates)) .* log(rand([1,1])); % exprnd_local(1/sum(rate_list));
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
            m(3) = m(3) + 1;
        elseif rxn == 4
            m(1) = m(1) - 1;
        elseif rxn == 5
            m(2) = m(2) - 1;
        elseif rxn == 6
            m(3) = m(3) - 1;
        elseif rxn == 7
            p(1) = p(1) + 1;
        elseif rxn == 8
            p(2) = p(2) + 1;
        elseif rxn == 9
            p(3) = p(3) + 1;
        elseif rxn == 10
            p(1) = p(1) - 1;
        elseif rxn == 11
            p(2) = p(2) - 1;
        elseif rxn == 12
            p(3) = p(3) - 1;
        elseif rxn == 13
            E = E-1;
            M = M+1;
        elseif rxn == 14
            E = E+1;
            M = M-1;
        end        
          
    else
        ref = step;
        break
        
    end
    
    alpha_x = 0.1 * (1 - p(end)/sum(p));
    
    t_vec_x(step) = t;
    EM_vec_x(step,:) = [E,M];
    m_vec_x(step,:) = m; 
    p_vec_x(step,:) = p;          
end

%%%%%
toc
%%%%%

%%% Post-processing

% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec_x(ref:end) = [];
    EM_vec_x(ref:end, :) = [];
    m_vec_x(ref:end,:) = [];
    p_vec_x(ref:end,:) = [];   
    rates_vec(ref:end,:) = [];
end

prot_frac_x = p_vec_x ./ sum(p_vec_x, 2);

%%%%%
toc
%%%%%



%% Sim MORE burden

%%% Starting variables
t = 0; % Time
E = N; % Engineered cell
M = E-N; % Mutant cell
m = [3,3,3]; % G,H
p = [1200,1200,1200]; % R,G,H

%%% Store vectors
t_vec = zeros(num_steps, 1); % Time vector
m_vec = zeros(num_steps, 3);
p_vec = zeros(num_steps, 3);
EM_vec = zeros(num_steps, 2); % E/M vector

%%% Record values at time 1
EM_vec(1,:) = [E,M];
p_vec(1,:) = p;
m_vec(1,:) = m;
alpha = 1 - p_vec(1,end) / sum(p_vec(1,:)); % Burden boost

%%%%%
tic
%%%%%

%%%
for step = 2:num_steps+1 % So we record first results as 0
       
%     step   

    % Assume all species equal / rates are molecs/h for a single
    % transcript - dynamcis/trends may vary if we scale to a whole cell
    
    %%% mRNA production
    gamma_km = (51-25*(E/(E+M))) / 26;
    km_R = km_tot/51 * gamma_km;
    km_G = km_tot*25/51 * gamma_km;
    km_H = km_tot*25/51 * (E/(E+M)); % Scale km by proportion of M!
    
    %%% mRNA decay   
    dm_R = m(1) * dm_tot/51;
    dm_G = m(2) * dm_tot*25/51;
    dm_H = m(3) * dm_tot*25/51;    
    
    %%% Protein production
    if p(1) > 0 && m(1)+m(2)+m(3) > 0
        gamma_kp = 51 / (p(1)*(m(1) + 25*m(2) + 25*m(3)));        
    else
        gamma_kp = 0;
    end
    kp_R = m(1) * p(1) * kp_tot/51 * gamma_kp;
    kp_G = m(2) * p(1) * kp_tot*25/51 * gamma_kp;
    kp_H = m(3) * p(1) * kp_tot*25/51 * gamma_kp;
    
    %%% Protein decay
    dp_R = p(1) * dp_tot/51;
    dp_G = p(2) * dp_tot*25/51;
    dp_H = p(3) * dp_tot*25/51;
    
    %%% Cell division
    lamEtoM = lam0 * (((E*z/(1+alpha)) + M) * E/(E+M)) ;%* (1+alpha); % E --> M
    lamMtoE = lam0 * (E*(1-z)) * M/(E+M) / (1+alpha); % M --> E
    
    % Record
    rates = [(E+M)*[km_R, km_G, km_H, dm_R, dm_G, dm_H, ...
                    kp_R, kp_G, kp_H, dp_R, dp_G, dp_H],...
                    lamEtoM, lamMtoE];                 
    rates_vec(step-1, :) = [rates(1:num_rxns-2)/(E+M), rates(num_rxns-1:num_rxns)];
    
    if E > 0 
        
        % Calculate next reaction and time from random numbers
        r1 = rand * sum(rates);
        r2 = -(1/sum(rates)) .* log(rand([1,1])); % exprnd_local(1/sum(rate_list));
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
            m(3) = m(3) + 1;
        elseif rxn == 4
            m(1) = m(1) - 1;
        elseif rxn == 5
            m(2) = m(2) - 1;
        elseif rxn == 6
            m(3) = m(3) - 1;
        elseif rxn == 7
            p(1) = p(1) + 1;
        elseif rxn == 8
            p(2) = p(2) + 1;
        elseif rxn == 9
            p(3) = p(3) + 1;
        elseif rxn == 10
            p(1) = p(1) - 1;
        elseif rxn == 11
            p(2) = p(2) - 1;
        elseif rxn == 12
            p(3) = p(3) - 1;
        elseif rxn == 13
            E = E-1;
            M = M+1;
        elseif rxn == 14
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

%%%%%
toc
%%%%%

%%% Post-processing

% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec(ref:end) = [];
    EM_vec(ref:end, :) = [];
    m_vec(ref:end,:) = [];
    p_vec(ref:end,:) = [];   
    rates_vec(ref:end,:) = [];
end

prot_frac = p_vec ./ sum(p_vec, 2);

%%%%%
toc
%%%%%



%% Plot LESS burden

% close all
% set(gcf, 'Position',  [0, 400, 800, 640])
% 
% %%% E/M
% plot(t_vec_x, EM_vec_x(:,1), '--', 'Color', E_yellow)
% hold on
% plot(t_vec_x, EM_vec_x(:,2), '--', 'Color', M_orange)
% % legend('E', 'M', 'Location', 'East')
% xlabel('Time / arb. unit')
% ylabel('Number of cells')
% grid on
% axis square
  
% plot(t_vec_x, prot_frac_x(:,1), 'Color', R_blue)
% hold on
% plot(t_vec_x, prot_frac_x(:,2), 'Color', G_green)
% plot(t_vec_x, prot_frac_x(:,3), 'Color', H_yellow)
% % legend('R','G','H', 'Location', 'East')
% xlabel('Time / arb. unit')
% ylabel('Protein fractions')
% grid on
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))



%% Plot MORE burden

close all
set(gcf, 'Position',  [0, 400, 800, 640])
% 
% %%% E/M
% plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
% hold on
% plot(t_vec, EM_vec(:,2), 'Color', M_orange)
% % legend('E', 'M', 'Location', 'East')
% xlabel('Time / arb. unit')
% ylabel('Number of cells')
% grid on
% axis square

%%% mRNA
plot(t_vec, m_vec(:,2), 'Color', G_green, 'LineWidth', 1.5)
hold on
plot(t_vec, m_vec(:,1), 'Color', R_blue, 'LineWidth', 1.5)
plot(t_vec, m_vec(:,3), 'Color', H_yellow, 'LineWidth', 1.5)
% legend('Genomic','Plasmid', 'Location', 'East')
xlabel('Time / arb. unit')
ylabel('Relative mRNA')
grid on
axis square

xlim([0 11.5])
ylim([0 20])
xticks([0, 2.5, 5.0, 7.5, 10])
set(gca,'Xticklabel',{'0','1','2', '3', '4'});

%%% Prot frac
% plot(t_vec, prot_frac(:,1), 'Color', R_blue, 'LineWidth', 1.5)
% hold on
% plot(t_vec, prot_frac(:,2), 'Color', G_green, 'LineWidth', 1.5)
% plot(t_vec, prot_frac(:,3), 'Color', H_yellow, 'LineWidth', 1.5)
% % legend('Ribosome','Genomic','Plasmid', 'Location', 'East')
% xlabel('Time / arb. unit')
% ylabel('Fraction of proteome')
% grid on
% axis square
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
% 
% xlim([0 11.5])
% % ylim([0 20])
% xticks([0, 2.5, 5.0, 7.5, 10])
% set(gca,'Xticklabel',{'0','1','2', '3', '4'});


%% Plot BOTH

close all
set(gcf, 'Position',  [0, 400, 800, 640])

%%% E/M
plot(t_vec_x, EM_vec_x(:,1), '--', 'Color', E_yellow)
hold on
plot(t_vec_x, EM_vec_x(:,2), '--', 'Color', M_orange)
plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
plot(t_vec, EM_vec(:,2), 'Color', M_orange)
% legend('E (low burden)', 'M (low burden)','E (high burden)','M (high burden)', 'Location', 'East')
% reorderLegend([2,1,4,3])
% xlim([0 8])
xlabel('Time / arb. unit')
ylabel('Number of cells')
grid on
axis square

%%% Prot frac
% plot(t_vec_x, prot_frac_x(:,1), 'Color', R_blue)
% hold on
% plot(t_vec_x, prot_frac_x(:,2), 'Color', G_green)
% plot(t_vec_x, prot_frac_x(:,3), 'Color', H_yellow)% 
% plot(t_vec, prot_frac(:,1), 'Color', R_blue)
% plot(t_vec, prot_frac(:,2), 'Color', G_green)
% plot(t_vec, prot_frac(:,3), 'Color', H_yellow)% 
% % legend('R (low burden)','G (low burden)','H (low burden)','R (high burden)','G (high burden)', 'H (high burden)','Location', 'East')
% % reorderLegend([])% 
% xlabel('Time / arb. unit')
% ylabel('Protein fractions')
% grid on
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))