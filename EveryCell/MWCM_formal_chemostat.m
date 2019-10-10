
% We know that alpha (and hence lam(M)) increases with more mutation. But
% by how much? This can be determined by running the WCM and getting a
% value for either growth rate or H fraction. For simplicity, here we are
% using the value of protein produced (which is a proxy for H fraction).

% NEW IDEA: instead of switching between models after every cell division
% (which varies), I could specify a time 't' to simulate divisions for,
% then pass over to my WCM dynamics for the current E/M values. By
% decreasing 't', I get closer to one division event happening, but it
% doesn't have to be 1 (could be 0 or 2+)?

% TO DO
% 2) Cells lose material too quickly - are rates correct?
% 6) Pre-allocate cell arrays and access last item properly...



%%

%%%%%
tic
%%%%%

close all; clear; clc;
set(0,'DefaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',18);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
c1 = [0.85 0.85 0.85];
c2 = [0.55 0.55 0.55];
c3 = [0.1 0.1 0.1];



%% Parameter set-up

%%% Global variables
N = 5e1; % Chemostat limit
steps = 2e8; % Steps in one simulation
% beta = 1; % Temporary

%%% Starting variables
t = 0; % Time
E0 = N; % Engineered cell / Total cells
M0 = N-E0; % Mutant cell
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate
m = [3,3]; % G,H
p = [3600,3600]; % G,H
num_species = 5; % t/m/m/p/p
n = 4.615; % Parameter in PDD equation
beta = 0.0001*exp(n*sum(p)/7200) - 0.0001; % Protein-based division (PBD) factor

%%% Cell population
pop = cell(1,N);

% Pre-allocate each cell with zeros - not yet fixed...
for i = 1:N
    pop{i} = zeros(5e5, num_species); % per cell rxn limit = 5e5, can change
end
cell_rxn_count = zeros(N,1); % Count how many rxns each cell has done for correct indexing

%%% Determine whether each cell is E/M
M_marker = zeros(N,1);

%%% Track time of E, M and time
t_vec = zeros(steps,1); % Could be useful
% EM_info = zeros(steps,3); % For recording only when a cell divides
EM_info = [0,E0,M0];

%%% Pre-populate ALL cells, plus E/M
for i = 1:N
    pop{i}(1,:) = [t, m, p];
end

% Variables to use in simulation
E = E0;
M = M0;

% No. possible reactions per cell
num_rxns = 10;

% Test - what is the average no. proteins when a cell divides?
% Would like it to be close to p-bar, hence I should adjust beta...
p_mean = [];
div_marker = zeros(N,1); % Check whether cell has already divided



%% Set up initial rate list

rate_list = zeros(1,num_rxns*N);

% mRNA and protein rates are in molecs/h. We assume that the rates
% for genomic and non-genomic gene expression are split evenly

km_G_0 = 24/2; km_H_0 = 24/2; % mRNA production = km = 24
dm_G_0 = 8/2; dm_H_0 = 8/2; % mRNA degradation = dm = 8
kp_G_0 = 24/2; kp_H_0 = 24/2; % Protein production = kp = 24
dp_G_0 = 0.02/2; dp_H_0 = 0.02/2; % Protein degradation = dp = 0.02
lamE_0 = lam0*(1-z);%*beta; % E production = lamE = 2*(1-z)*beta
lamM_0 = lam0*z*beta; % M production = lamM = 2*z*beta

% If M: km=0, lamE=0
    
% Iterate over each cell to complete the rate list
for k = 1:N    
    rate_list(num_rxns*(k-1)+1: k*num_rxns) = ...
        [km_G_0, km_H_0, dm_G_0, dm_H_0, kp_G_0, kp_H_0, dp_G_0, dp_H_0, lamE_0, lamM_0];
end



%% Gillespie steps

for i = 2:steps+1 % So we record first results as 0
    
%     i                          
   
    %%% CHOOSE REACTION %%%
    if E > 0 % Mutation can spread more
        
        % Get random numbers
        r1 = rand * sum(rate_list);
        r2 = -(1/sum(rate_list)) .* log(rand([1,1])); % exprnd_local(1/sum(rate_list));
        t = t + r2; % Update time
        
        % Get index of next reaction in rate_list
        % r1 < cumsum(rate_list): 0 when false, 1 when true
        % sum(x==0): count the zeros... +1 for first rxn
        rxn_index = sum(r1 < cumsum(rate_list)==0) + 1;
        
        % The cell no. (/N)
        mycell = ceil(rxn_index/num_rxns);
        cell_rxns = cell_rxn_count(mycell); % For indexing
        
        % The reaction no. (/10)
        rxn = mod(rxn_index, num_rxns);
        if rxn == 0
            rxn = num_rxns; % Convert 0 to 10 to help intuition
        end
                       
        % Previous values for our cell, for ease
        m = pop{mycell}(cell_rxns+1,2:3);
        p = pop{mycell}(cell_rxns+1,4:5);
        
        %%% UPDATE VALUES %%%
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
            
            % N+1th cell            
            E = E+1;
            m_div = floor(m/2); % Daughter loses half m/p (rounded down)
            p_div = floor(p/2);
            
            % Kick out cell randomly
            ind = randi(E+M);            
            
            % If it's the daughter, we just lose E
            if ind == E+M 
                E = E-1;       
                
            % Otherwise, our daughter replaces the one we lost
            else 
                % Update E/M
                if M_marker(ind) == 0
                    E = E-1;
                else
                    M = M-1;
                end
                
                % No. rxns of daughter (ind)
                cell_rxns_div = cell_rxn_count(ind);
                
                % New E replaces ind
                % Cells start with 1 row (+1) / we're recording for the next rxn (+1)
                pop{ind}(cell_rxns_div+2,:) = [t, m_div, p_div];
                
                % Record this cell as E
                M_marker(ind) = 0;
                
                % Update rate_list using daughter's rxn indices
                cellstart = num_rxns*(ind-1) + 1;
                cellend = num_rxns*(ind-1) + num_rxns;
                beta = 0.0001*exp(n*sum(p_div)/7200) - 0.0001; % PBD factor
                
                rate_list(cellstart : cellend) = [24/2, 24/2, m_div(1)*8/2, m_div(2)*8/2,...
                                                  m_div(1)*24/2, m_div(2)*24/2, p_div(1)*0.02/2, p_div(2)*0.02/2,...
                                                  lam0*(1-z)*beta, lam0*z*beta];
            end         
            
            % Parent cell
            % If we lose this one (ind = mycell), no separate split in mass is needed
            
            %%%%% TESTING
            % Check that cell has divided at least x* previously
            if div_marker(mycell) >= 1
                p_mean(end+1,:) = sum(p); % Record protein level upon division
            end
            %%%%%
            
            if ind ~= mycell
                m = ceil(m/2); % Parent loses half m/p (rounded up)
                p = ceil(p/2);
            end
            
            % Record E/M with time
            EM_info(end+1,:) = [t,E,M];
            
            %%%%% TESTING
            div_marker(mycell) = div_marker(mycell) + 1;
            %%%%%
        
        elseif rxn == 10
            
            % N+1th cell            
            M = M+1;
            m_div = floor(m/2); % Daughter loses half m/p (rounded down)
            p_div = floor(p/2);
            
            % Kick out cell randomly
            ind = randi(E+M);
            
            % If it's the daughter, we just lose M
            if ind == E+M 
                M = M-1;       
                
            % Otherwise, our daughter replaces the one we lost
            else 
                % Update E/M
                if M_marker(ind) == 0
                    E = E-1;
                else
                    M = M-1;
                end
                
                % No. rxns of daughter (ind)
                cell_rxns_div = cell_rxn_count(ind);
                
                % New M replaces ind
                % Cells start with 1 row (+1) / we're recording for the next rxn (+1)
                pop{ind}(cell_rxns_div+2,:) = [t, m_div, p_div];
                
                % Record this cell as M
                M_marker(ind) = 1;
                
                % Update rate_list using daughter's rxn indices
                cellstart = num_rxns*(ind-1) + 1;
                cellend = num_rxns*(ind-1) + num_rxns;
                beta = 0.0001*exp(n*sum(p_div)/7200) - 0.0001; % PBD factor
                
                rate_list(cellstart : cellend) = [24/2, 0, m_div(1)*8/2, m_div(2)*8/2,...
                                                  m_div(1)*24/2, m_div(2)*24/2, p_div(1)*0.02/2, p_div(2)*0.02/2,...
                                                  lam0*0*beta, lam0*1*beta];
            end         
            
            % Parent cell
            % If we lose this one (ind = mycell), no separate split in mass is needed
            
            %%%%% TESTING
            % Check that cell has divided at least x* previously
            if div_marker(mycell) >= 1
                p_mean(end+1,:) = sum(p); % Record protein level upon division
            end
            %%%%%
            
            if ind ~= mycell
                m = ceil(m/2); % Parent loses half m/p (rounded up)
                p = ceil(p/2);
            end
            
            % Record E/M with time
            EM_info(end+1,:) = [t,E,M];
            
            %%%%% TESTING
            div_marker(mycell) = div_marker(mycell) + 1;
            %%%%%
        end            

    else % If sim step lim reached
        ref = i;
        break
    end

    %%% RECORD FOR AFFECTED CELL %%%    
    % Store time
    t_vec(i) = t;   
    
    % Species numbers
    % Cells start with 1 row (+1) / we're recording for the next rxn (+1)
    pop{mycell}(cell_rxns+2,:) = [t, m, p];
    
    % Update no. rxns for our cell (for indexing)
    cell_rxn_count(mycell) = cell_rxn_count(mycell) + 1;
    
    % Update rate_list using parent's rxn indices 
    cellstart = rxn_index - rxn + 1;
    cellend = rxn_index - rxn + num_rxns;
    beta = 0.0001*exp(n*sum(p)/7200) - 0.0001; % PBD factor
    
    if M_marker(mycell) == 0 % Cell = E
        rate_list(cellstart : cellend) = [24/2, 24/2, m(1)*8/2, m(2)*8/2,...
                                          m(1)*24/2, m(2)*24/2, p(1)*0.02/2, p(2)*0.02/2,...
                                          lam0*(1-z)*beta, lam0*z*beta];            
    else % Cell = M
        rate_list(cellstart : cellend) = [24, 0, m(1)*8/2, m(2)*8/2,...
                                          m(1)*24/2, m(2)*24/2, p(1)*0.02/2, p(2)*0.02/2,...
                                          lam0*0*beta, lam0*1*beta];
    end
      
end

% If we reached the end early, trim the vectors!
if i < steps+1
    t_vec(ref:end) = [];
    EM_info(ref:end, :) = [];
end

%%%%%
toc
%%%%%

% Remove trailing zeros from each cell
for q = 1:N
    num = cell_rxn_count(q)+2; % Index of the first 0
    pop{q}(num:end, :) = []; % Delete all these
end

disp(mean(p_mean));
disp(std(p_mean));



%% E/M over time

set(gcf, 'Position',  [1200, 4, 800, 640])

plot(EM_info(:,1), EM_info(:,2), 'Color', E_yellow)
hold on
plot(EM_info(:,1), EM_info(:,3), 'Color', M_orange)
l = legend('E','M','Location','East');
xlabel('Time / h')
ylabel('Cell count')
% ylim([0 1e4])
grid on

% set(gca, 'Color', 'none');
% set(l,'color','none');
% export_fig MWCMformal_poster1.pdf -transparent



%% Old code

% Get index of reaction
% for j = 1:length(rate_list)
%     % As soon as r1 doesn't pass the cumulative sum, that's our rxn
%     if r1 < sum(rate_list(1:j))
%         rxn_index = j;
%         break
%     end
% end
