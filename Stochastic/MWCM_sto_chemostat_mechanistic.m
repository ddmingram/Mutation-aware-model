
%%% To Do

% Will probs need to change km_H decay... make more natural


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
N = 1e2; % Chemostat limit
num_steps = 5e8; % Steps in one simulation
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate
num_rxns = 20; % No. possible reactions per cell

%%% Starting variables
% Core
t = 0; % Time
E = N; % Engineered cell
M = E-N; % Mutant cell
% Transcription
D_R = 100; D_G = 100; D_H = 100;
pol = 1000;
TX_R = 0; TX_G = 0; TX_H = 0;
m = [3,3,3];
% Translation
TL_R = 0; TL_G = 0; TL_H = 0;
p = [100,1.4e5,1.4e5];

%%% Store vectors
t_vec = zeros(num_steps, 1); % Time vector
m_vec = zeros(num_steps, 3);
p_vec = zeros(num_steps, 3);
EM_vec = zeros(num_steps, 2); % E/M vector

%%% Record values at time 1
EM_vec(1,:) = [E,M];
p_vec(1,:) = p;
m_vec(1,:) = m;

%%% Burden
A = 1; % Burden effect
alpha = A * (p_vec(1,end) / sum(p_vec(1,:)));

%%% Relative base rates
kTX = 1;
km = 1;
dm = 1/3;
kTL = 1;
kp = 1;
dp = 1/1200;

kp_tot = 240;

%%% Test
% track_rates = zeros(num_steps, num_rxns);
% track_species = zeros(num_steps, 16);

%%%%%
tic
%%%%%

%%
for step = 2:num_steps+1 % So we record first results as 0
       
%     step   

    % Ribosomal synthesis is 1/25th of average
    % Relative rates based on molecs/h for a single transcript 
    % Assume no decay of D/pol/TX/TL
    
    %%%%%%% Transcription %%%%%%%
    
    %%% DNA-pol complex production
    kTX_R = D_R * pol * kTX;
    kTX_G = D_G * pol * kTX;
    kTX_H = D_H * pol * kTX * (E/(E+M));  
    
    %%% mRNA production
    km_R = TX_R * km/25;
    km_G = TX_G * km;
    km_H = TX_H * km;
    
    %%% Decay  
    dm_R = m(1) * dm/25;
    dm_G = m(2) * dm;
    dm_H = m(3) * dm;
    
    %%%%%%% Translation %%%%%%%
    
    %%% mRNA-ribosome complex production
    kTL_R = m(1) * p(1) * kTL;
    kTL_G = m(2) * p(1) * kTL;
    kTL_H = m(3) * p(1) * kTL;    
    
    %%% Protein production    
%     kp_R = TL_R * kp/25;
%     kp_G = TL_G * kp;
%     kp_H = TL_H * kp;    
    if TL_R + TL_G + TL_H > 0
        gamma = kp_tot / (kp * (TL_R/25 + TL_G + TL_H));
    else
        gamma = 0;
    end
    kp_R = TL_R * kp/25 * gamma;
    kp_G = TL_G * kp * gamma;
    kp_H = TL_H * kp * gamma;
    
    %%% Decay
    dp_R = p(1) * dp/25;
    dp_G = p(2) * dp;
    dp_H = p(3) * dp;
    
    %%%%%%% Cell division %%%%%%%
    
    %%% Cell division
    lamEtoM = lam0 * (E*z/(1+alpha) + M) * E/(E+M); % E --> M
    lamMtoE = lam0 * (E*(1-z)/(1+alpha)) * M/(E+M); % M --> E
    
    % Record
    rates = [ N*[kTX_R, kTX_G, kTX_H, km_R, km_G, km_H, dm_R, dm_G, dm_H, ...
                 kTL_R, kTL_G, kTL_H, kp_R, kp_G, kp_H, dp_R, dp_G, dp_H],...
                 lamEtoM, lamMtoE];
    
%     track_rates(step-1,:) = rates;    
%     track_species(step-1,:) = [D_R, D_G, D_H, pol,...
%                                TX_R, TX_G, TX_H, m,...
%                                TL_R, TL_G, TL_H, p];
                       
    if E > 0 
        
        % Calculate next reaction and time from random numbers
        r1 = rand * sum(rates);
        r2 = -(1/sum(rates)) .* log(rand([1,1])); % exprnd(1/sum(rate_list));
        t = t + r2;
        
        % Get index of next reaction in rate_list
        % r1 < cumsum(rates): 0 when false, 1 when true
        % sum(x==0): count the zeros... +1 for first rxn
        rxn = sum(r1 < cumsum(rates)==0) + 1;     
        
        % TX complex form
        if rxn == 1
            D_R = D_R - 1;
            pol = pol - 1;
            TX_R = TX_R + 1;            
        elseif rxn == 2
            D_G = D_G - 1;
            pol = pol - 1;
            TX_G = TX_G + 1;            
        elseif rxn == 3
            D_H = D_H - 1;
            pol = pol - 1;
            TX_H = TX_H + 1;
        
        % mRNA form
        elseif rxn == 4
            TX_R = TX_R - 1;
            D_R = D_R + 1;
            pol = pol + 1;
            m(1) = m(1) + 1;            
        elseif rxn == 5
            TX_G = TX_G - 1;
            D_G = D_G + 1;
            pol = pol + 1;
            m(2) = m(2) + 1;       
        elseif rxn == 6
            TX_H = TX_H - 1;
            D_H = D_H + 1;
            pol = pol + 1;
            m(3) = m(3) + 1;
        
        % mRNA decay
        elseif rxn == 7
            m(1) = m(1) - 1;        
        elseif rxn == 8
            m(2) = m(2) - 1;         
        elseif rxn == 9
            m(3) = m(3) - 1;
        
        % TL complex form
        elseif rxn == 10
            m(1) = m(1) - 1;
            p(1) = p(1) - 1;
            TL_R = TL_R + 1;        
        elseif rxn == 11
            m(2) = m(2) - 1;
            p(1) = p(1) - 1;
            TL_G = TL_G + 1;        
        elseif rxn == 12
            m(3) = m(3) - 1;
            p(1) = p(1) - 1;
            TL_H = TL_H + 1;
        
        % Protein form
        elseif rxn == 13
            TL_R = TL_R - 1;
            m(1) = m(1) + 1;
            p(1) = p(1) + 2;        
        elseif rxn == 14
            TL_G = TL_G - 1;
            m(2) = m(2) + 1;
            p(1) = p(1) + 1;
            p(2) = p(2) + 1;        
        elseif rxn == 15
            TL_H = TL_H - 1;
            m(3) = m(3) + 1;
            p(1) = p(1) + 1;
            p(3) = p(3) + 1;
        
        % Protein decay
        elseif rxn == 16
            p(1) = p(1) - 1;       
        elseif rxn == 17
            p(2) = p(2) - 1;        
        elseif rxn == 18
            p(3) = p(3) - 1;
        
        % Cell division
        elseif rxn == 19
            E = E-1;
            M = M+1;        
        elseif rxn == 20
            E = E+1;
            M = M-1;
        end        
          
    else
        ref = step;
        break
        
    end
    
    alpha = A * (p(end)/sum(p));
    
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
%     track_rates(ref:end,:) = [];
%     track_species(ref:end,:) = [];
end

% How alpha caries over time
alpha_vec = zeros(size(p_vec,1), 1);
for i = 1:length(alpha_vec)
    alpha_vec(i) = A * (p_vec(i,end)/sum(p_vec(i,:)));
end


%%%%%
toc
%%%%%


%% EM plot

close all
set(gcf, 'Position',  [1600, 400, 800, 640])

plot(t_vec, EM_vec(:,1), 'Color', E_yellow)
hold on
plot(t_vec, EM_vec(:,2), 'Color', M_orange)
legend('E', 'M')
xlabel('Time / h')
ylabel('Number of cells')
% xlim([0 8])
grid on
axis square



%% Plot mRNA

% close all
% set(gcf, 'Position',  [1600, 400, 800, 640])
% 
% plot(t_vec, m_vec(:,1), 'LineWidth', 1.5, 'Color', R_blue)
% hold on
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

% close all
% set(gcf, 'Position',  [1600, 400, 800, 640])

