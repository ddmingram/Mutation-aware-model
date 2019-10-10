
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
R_blue = [0 114 189]/255;
G_green = [119 172 48]/255;
H_yellow = [237 177 32]/255;


%% Variables

%%% Global variables
N = 1e4; % Chemostat limit
num_steps = 1e8; % Steps in one simulation
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate

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

%%% Useful other
num_rxns = 14; % No. possible reactions per cell

%%% NEW TESTING
km_tot = 24;
dm_tot = 8;
kp_tot = 24;
dp_tot = 0.02;
rates_vec = zeros(num_steps,14);

%%%%%
tic
%%%%%

%%
for step = 2:num_steps+1 % So we record first results as 0
       
%     step   

    % Assume all species equal / rates are molecs/h for a single
    % transcript - dynamcis/trends may vary if we scale to a whole cell
    
    %%% mRNA production
    gamma_km = (51-25*(E/(E+M))) / 26;
    km_R = km_tot/51 * gamma_km;
    km_G = km_tot*25/51 * gamma_km;
    km_H = km_tot*25/51 * (E/(E+M)); % Scale km by proportion of M!

%     km_G = 24/2; %%% No scale
%     km_H = 24/2 * (E/(E+M));

%     gamma_km = 2 / ((1 + E/(E+M))); %%% Scale all
%     km_G = 12 * gamma_km;
%     km_H = 12 * E/(E+M) * gamma_km;
    
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
    
%     kp_R = m(1) * p(1) * kp_tot/51;
%     kp_G = m(2) * p(1) * kp_tot*25/51;
%     kp_H = m(3) * p(1) * kp_tot*25/51;
    
    %%% Protein decay
    dp_R = p(1) * dp_tot/51;
    dp_G = p(2) * dp_tot*25/51;
    dp_H = p(3) * dp_tot*25/51;
    
    %%% Cell division
    lamEtoM = lam0 * ((E*z+M) * E/(E+M)) * (1+alpha); % E --> M
    lamMtoE = lam0 * (E*(1-z)) * M/(E+M); % M --> E
    
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

% test_vec = [rates_vec(:,7), rates_vec(:,8)];
% plot(sum(test_vec,2))


%% Post-processing

% If we reached the end early, trim the vectors!
if step < num_steps+1
    t_vec(ref:end) = [];
    EM_vec(ref:end, :) = [];
    m_vec(ref:end,:) = [];
    p_vec(ref:end,:) = [];   
    rates_vec(ref:end,:) = [];
end

% How alpha caries over time
alpha_vec = zeros(size(p_vec,1), 1);
for i = 1:length(alpha_vec)
    alpha_vec(i) = 1 - ( p_vec(i,end)/sum(p_vec(i,:)) );
end

% For plotting transition rate
lamEtoM = lam0 * ((EM_vec(:,1).*z + EM_vec(:,2)) .* EM_vec(:,1)./(EM_vec(:,1)+EM_vec(:,2))).*(1+alpha_vec);


%%%%%
toc
%%%%%


%% EM plot

% set(gcf, 'Position',  [0, 400, 800, 640])
% 
% plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
% hold on
% plot(t_vec, EM_vec(:,2), 'Color', M_orange)
% legend('E', 'M')
% xlabel('Time / h')
% ylabel('Number of cells')
% % xlim([0 8])
% grid on



%% EM plot + H fraction

close all
set(gcf, 'Position',  [0, 400, 800, 640])

% subplot(1,2,1)

plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
hold on
plot(t_vec, EM_vec(:,2), 'Color', M_orange)
legend('E', 'M')
xlabel('Time / arb. unit')
ylabel('Number of cells')
% xlim([0 8])
grid on
axis square

% subplot(1,2,2)
% 
% plot(t_vec, 1-alpha_vec, 'LineWidth', 0.5, 'Color',mediumgrey)
% xlabel('Time / h')
% ylabel('H fraction')
% % xlim([0 8])
% grid on
% axis square
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))



%% Plot mRNA

% close all
% set(gcf, 'Position',  [1600, 400, 800, 640])
% 
% % plot(t_vec, m_vec(:,1), 'LineWidth', 1.5, 'Color', R_blue)
% 
% plot(t_vec, m_vec(:,2), 'LineWidth', 1.5, 'Color', G_green)
% hold on
% plot(t_vec, m_vec(:,3), 'LineWidth', 1.5, 'Color', H_yellow)
% xlabel('Time / h')
% ylabel('Quantity of mRNA')
% legend('G','H', 'Location', 'East')
% grid on
% axis square


%% Plot proteins

% close all
% set(gcf, 'Position',  [1600, 400, 800, 640])
% 
% plot(t_vec(1:step-1), p_vec(1:step-1,1), 'LineWidth', 1.5, 'Color', R_blue)
% hold on
% plot(t_vec(1:step-1), p_vec(1:step-1,2), 'LineWidth', 1.5, 'Color', G_green)
% plot(t_vec(1:step-1), p_vec(1:step-1,3), 'LineWidth', 1.5, 'Color', H_yellow)
% xlabel('Time / h')
% ylabel('Quantity of protein')
% legend('R','G','H', 'Location', 'East')
% % xlim([0 8])
% % ylim([0 2e6])
% grid on
% axis square


%% Other testing

% set(gcf, 'Position',  [1600, 400, 800, 640])
% 
% plot(t_vec, rates_vec(:,3))
% hold on
% plot(t_vec, rates_vec(:,4))
% grid on
% axis square
