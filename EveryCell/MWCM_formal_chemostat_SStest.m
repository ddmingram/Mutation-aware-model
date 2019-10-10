
% We know that alpha (and hence lam(M)) increases with more mutation. But
% by how much? This can be determined by running the WCM and getting a
% value for either growth rate or H fraction. For simplicity, here we are
% using the value of protein produced (which is a proxy for H fraction).

% TO DO
% 1) Standard simulations go above cell max and below 0 at the end??
% 2) Cells lose material too quickly - are rates correct?
% 3) For recording other variables, can either record straight after they
% change (like with E/M), or at end of each Gillespie step for each cell
% 4) Something not right with cell replacement - I end up with cells but
% with nothing in them for majority of time steps?
% 5) REF always encountered at last step no matter what??

%%% Steady state test:
% Just E cells, no division, mRNA/proteins start at 0. See what levels the
% species reach naturally.


%%

%%%%%
tic
%%%%%

close all; clear; clc;
set(0,'DefaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',18);

c1 = [0.85 0.85 0.85];
c2 = [0.55 0.55 0.55];
c3 = [0.1 0.1 0.1];


%%

%%% Global variables
N = 10; % Chemostat limit
steps = 2e6; % Steps in one simulation

%%% Starting variables
t = 0; % Time
E0 = N; % Engineered cell / Total cells
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate
m = [0,0]; % G,H
p = [0,0]; % G,H

%%% Cell population
pop = cell(1,N);

%%% Track time of E, M and time
t_vec = zeros(steps,1); % Could be useful

%%% Pre-populate ALL cells, plus E/M
for i = 1:N
    pop{i}(1,:) = [t, m, p];
end

% Variables to use in simulation
E = E0;

% No. possible reactions per cell
rxns = 8;


%%

for i = 2:steps+1 % So we record first results as 0
    
%     i

    rate_list = zeros(1,rxns*N);
    
    % ITERATE OVER EACH CELL TO GET TOTAL RATES
    for k = 1:N
 
        % Access m/p for each cell, = previous step
        m = pop{k}(end,2:3);
        p = pop{k}(end,4:5);

        % Cell is either E or M, via M_marker=0/1

      
        % Assume that rates for genomic and non-genomic gene expression are
        % split evenly

        % mRNA production
        km_G = 24/2; % molecs/h (assuming all species equal)
        km_H = 24/2;

        % mRNA decay
        dm_G = m(1) * 8/2; % molecs/h (assuming all species equal)
        dm_H = m(2) * 8/2;    

        % Protein production
        kp_G = m(1) * 24/2; % molecs/h (assuming all species equal)
        kp_H = m(2) * 24/2;

        % Protein decay
        dp_G = p(1) * 0.02/2; % molecs/h (assuming all species equal)
        dp_H = p(2) * 0.02/2;  
        
        rate_list(8*(k-1)+1: k*8) = [km_G, km_H, dm_G, dm_H, kp_G, kp_H, dp_G, dp_H];                            
    end
    
    % CHOOSE REACTION AND UPDATE
    if E > 0 % Mutation can spread more
        
        % Get random numbers
        r1 = rand * sum(rate_list);
        r2 = exprnd(1/sum(rate_list));
        t = t + r2; % Update the time
        
        % Get index of next reaction in rate_list
        % r1 < cumsum(rate_list): 0 when false, 1 when true
        % sum(x==0): count the zeros... +1 for first rxn
        rxn_index = sum(r1 < cumsum(rate_list)==0) + 1;
        
        % The reaction no. (/10)
        rxn = mod(rxn_index, 8);
        % The cell no. (/N)
        mycell = ceil(rxn_index/8);
        
        % Previous values for our cell, for ease
        m = pop{mycell}(end,2:3);
        p = pop{mycell}(end,4:5);
        
        % Update values
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
        elseif rxn == 0
            p(2) = p(2) - 1;            
        end            

    else % If sim step lim reached
        ref = i;
        break
    end

    %%% Update info!

    t_vec(i) = t; % Useful
    
    pop{mycell}(end+1,:) = [t, m, p];

end


%%%%%
toc
%%%%%

% For 2e6 steps, stationary distribution starts from ~60,000. Start from
% 70,000 for good measure.

GmRNAmeans = zeros(1,N);
HmRNAmeans = zeros(1,N);
Gproteinmeans = zeros(1,N);
Hproteinmeans = zeros(1,N);

for i = 1:N
    GmRNAmeans(i) = mean(pop{i}(7e4:end,2));
    HmRNAmeans(i) = mean(pop{i}(7e4:end,3));
    Hproteinmeans(i) = mean(pop{i}(7e4:end,4));
    Gproteinmeans(i) = mean(pop{i}(7e4:end,5));
end

disp( mean(GmRNAmeans) ) % G mRNA mean
disp( mean(HmRNAmeans) ) % H mRNA mean
disp( mean(Gproteinmeans) ) % G protein mean
disp( mean(Hproteinmeans) ) % H protein mean


%% E/M over time

set(gcf, 'Position',  [1200, 400, 800, 640])

plot(pop{2}(:,1), pop{2}(:,2))
hold on
plot(pop{2}(:,1), pop{2}(:,3))
l = legend('Genomic','Heterologous','Location','NorthEast');
xlabel('Time / h')
ylabel('Species count')
grid on
hold off

% To plot the mean graph, take the lowest number of steps from all the
% cells... ~190,000





