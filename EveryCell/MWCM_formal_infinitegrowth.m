
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
% 1) Trim BigArray to remove zeros before new cells start, OR replace with NAs
% 2) Cells lose material too quickly - are rates correct?


%%
close all; clear; clc;
set(0,'DefaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',18);

E_yellow = [0.88 0.78 0.02];
M_orange = [0.88 0.53 0];
c1 = [0.85 0.85 0.85];
c2 = [0.55 0.55 0.55];
c3 = [0.1 0.1 0.1];


%%

%%% Global variables
steps = 1e3; % Steps in one simulation
cells = 75000; % Max cells being considered (up to 75,000)

%%% Starting variables
t = 0; % Time
E = 1; % Engineered cell
M = 0; % Mutant cell
lam0 = 2; % Divisions / h
z = 0.01; % Mutation rate
m = [1,1,1,1]; % R,E,Q,H
p = [400,400,400,400]; % R,E,Q,H

%%% Big Store Array
% Rows = simulation steps
% Columns = m(4), p(4)
% Pages = cells
BigArray = zeros(steps, 9, cells);

%%% Determine whether each cell is E/M
M_marker = zeros(cells,1);

%%% Track time for E or M increases
EM_info = zeros(cells,2); % Time updates for cell x when any cell divides

%%% Populate info for first cell
BigArray(1,1,1) = t;
BigArray(1,2:5,1) = m;
BigArray(1,6:9,1) = p;



%%

for i = 2:steps+1 % So we record first results as 0
    
%     i
    
    % Stop after cell limit - INFINITEGROWTH
    if (E+M) < cells 
    
        % For each Gillespie step, iterate over every cell!
        for k = 1:(E+M)
            
%             k
            
            % Cell is either E or M, via M_marker=0/1

            % Access t/m/p for each cell, = previous step
            t = BigArray(i-1,1,k);
            m = BigArray(i-1,2:5,k);
            p = BigArray(i-1,6:9,k);            

            % Cell division PER CELL
            if M_marker(k) == 0 % If E, normal division
                lamE = lam0 * (1-z);
                lamM = lam0 * z;
            else % If M, any division is M
                lamE = lam0 * 0;
                lamM = lam0 * 1;
            end

            % mRNA production
            km_R = 24/4; % molecs/h (assuming all species equal)
            km_E = 24/4;
            km_Q = 24/4;
            if M_marker(k) == 0 % If E, normal H transcription
                km_H = 24/4;
            else % If M, no H transcription
                km_H = 0;
            end

            % mRNA decay
            dm_R = m(1) * 8/4; % molecs/h (assuming all species equal)
            dm_E = m(2) * 8/4;
            dm_Q = m(3) * 8/4;
            dm_H = m(4) * 8/4;    

            % Protein production
            kp_R = m(1) * 24/4; % molecs/h (assuming all species equal)
            kp_E = m(2) * 24/4;
            kp_Q = m(3) * 24/4;
            kp_H = m(4) * 24/4;

            % Protein decay
            dp_R = p(1) * 0.02/4; % molecs/h (assuming all species equal)
            dp_E = p(2) * 0.02/4;
            dp_Q = p(3) * 0.02/4;
            dp_H = p(4) * 0.02/4;  

            % PER CELL rates
            rates = [km_R, km_E, km_Q, km_H, dm_R, dm_E, dm_Q, dm_H,...
                     kp_R, kp_E, kp_Q, kp_H, dp_R, dp_E, dp_Q, dp_H,...
                     lamE, lamM];            

            % Calculate next reaction and time from random numbers
            r1 = rand * sum(rates);
            r2 = exprnd(1/sum(rates));
            t = t + r2;        

            if r1 < rates(1)
                m(1) = m(1) + 1;
            elseif r1 >= rates(1) && r1 < sum(rates(1:2))
                m(2) = m(2) + 1;
            elseif r1 >= sum(rates(1:2)) && r1 < sum(rates(1:3))
                m(3) = m(3) + 1;
            elseif r1 >= sum(rates(1:3)) && r1 < sum(rates(1:4))
                m(4) = m(4) + 1;
            elseif r1 >= sum(rates(1:4)) && r1 < sum(rates(1:5))
                p(1) = p(1) + 1;
            elseif r1 >= sum(rates(1:5)) && r1 < sum(rates(1:6))
                p(2) = p(2) + 1;
            elseif r1 >= sum(rates(1:6)) && r1 < sum(rates(1:7))
                p(3) = p(3) + 1;
            elseif r1 >= sum(rates(1:7)) && r1 < sum(rates(1:8))
                p(4) = p(4) + 1;
            elseif r1 >= sum(rates(1:8)) && r1 < sum(rates(1:9))
                m(1) = m(1) - 1;
            elseif r1 >= sum(rates(1:9)) && r1 < sum(rates(1:10))
                m(2) = m(2) - 1;
            elseif r1 >= sum(rates(1:10)) && r1 < sum(rates(1:11))
                m(3) = m(3) - 1;
            elseif r1 >= sum(rates(1:11)) && r1 < sum(rates(1:12))
                m(4) = m(4) - 1;
            elseif r1 >= sum(rates(1:12)) && r1 < sum(rates(1:13))
                p(1) = p(1) - 1;
            elseif r1 >= sum(rates(1:13)) && r1 < sum(rates(1:14))
                p(2) = p(2) - 1;
            elseif r1 >= sum(rates(1:14)) && r1 < sum(rates(1:15))
                p(3) = p(3) - 1;
            elseif r1 >= sum(rates(1:15)) && r1 < sum(rates(1:16))
                p(4) = p(4) - 1;

            elseif r1 >= sum(rates(1:16)) && r1 < sum(rates(1:17))
                E = E+1;
                % New cell gains half (rounded down). Start t = current t
                BigArray(i,:,(E+M)) = [t, floor([m,p]/2)];
                % Get time of E+1              
                EM_info(E+M,:) = [t,0]; 
                
            else
                M = M+1;
                % New cell gains half (rounded down). Start t = current t
                BigArray(i,:,(E+M)) = [t, floor([m,p]/2)]; % " "
                % Get time of M+1            
                EM_info(E+M,:) = [t,1];
                % 0-->1 = we have M
                M_marker(E+M) = 1;
            end

            % For each cell, update new values from Gillespie step
            if r1 < sum(rates(1:16))
                BigArray(i,:,k) = [t,m,p];      
            else % If division, then current cell loses half mRNA/protein
                BigArray(i,:,k) = [t, ceil([m,p]/2)];
            end
        end
        
    else % If cell lim reached
        ref = i;
        break
    end
end

% If we reached the end early, trim the vectors!
if i < steps+1
    BigArray(ref:end, :, :) = [];
end


%% Calculate EM increases vs. time using COARSE and FINE methods

%%% Time updates for cell x when any cell divides
EM_info = sortrows(EM_info); % Sort by first column

E_info = zeros(length(EM_info),2); % E-specific increase
E_info(:,1) = EM_info(:,1); % 1st col = time increases
M_info = zeros(length(EM_info),2); % " "
M_info(:,1) = EM_info(:,1);

E_count = 0; % Track whether a time point is E+1 or M+1, using EM_info_fine
M_count = 0;
for i = 1:length(EM_info)
    if EM_info(i,2) == 0 % If E cell
        E_count = E_count + 1;
        E_info(i,2) = E_count;
        M_info(i,2) = M_count;
    else % If M cell
        M_count = M_count + 1;
        E_info(i,2) = E_count;
        M_info(i,2) = M_count;
    end
end


%% EM plots

% set(gcf, 'Position',  [1200, 400, 800, 640])
% 
% plot(E_info(:,1), E_info(:,2), 'Color', E_yellow)
% hold on
% plot(M_info(:,1), M_info(:,2), 'Color', M_orange)
% xlabel('Time / h')
% ylabel('Number')
% legend('E', 'M', 'Location', 'North')
% 
% grid on



%% Specific cell plot - Protein

% set(gcf, 'Position',  [1200, 400, 800, 640])
% 
% plot(BigArray(:,1,1), BigArray(:,9,1))
% hold on
% plot(BigArray(266:end,1,49), BigArray(266:end,9,49)) % Change time points and cell no.
% xlabel('Time / h')
% ylabel('Quantity of H protein')
% legend('First E cell', 'First M cell', 'Location', 'East')
% grid on



%% Specific cell plot - mRNA

% set(gcf, 'Position',  [1200, 400, 800, 640])
% 
% plot(BigArray(:,1,1), BigArray(:,5,1))
% hold on
% plot(BigArray(266:end,1,49), BigArray(266:end,5,49)) % Change time points and cell no.
% xlabel('Time / h')
% ylabel('Quantity of mRNA_H')
% ylim([0 15])
% legend('First E cell', 'First M cell', 'Location', 'East')
% grid on








