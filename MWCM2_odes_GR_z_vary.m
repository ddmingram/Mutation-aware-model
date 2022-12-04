%% Care

% Any parameter vectors I use in ODEs must be column vectors, not row
% vectors, e.g. prom_vec and RBS_vec. E.g. prom_vec*plasmid*P =
% [a,b,c]*1*[x;y;z] = single value (3*3 matrix multplication).

function [dydt,y] = MWCM2_odes_GR_z_vary(t, y, params)

%% Unpack params, populate variables

%%%%%%%%%% PARAMETERS

x = 1;

n = params(x); x=x+1;
N = params(x); x=x+1;
mass = params(x); x=x+1;
nut = params(x); x=x+1;
nq = params(x); x=x+1;
v_e = params(x); x=x+1;
K_e = params(x); x=x+1;
n_R = params(x); x=x+1;
n_C = params(x); x=x+1;
n_Q = params(x); x=x+1;
n_H = params(x); x=x+1;
v_TX_R = params(x); x=x+1;
v_TX_C = params(x); x=x+1;
v_TX_Q = params(x); x=x+1;
K_TX_R = params(x); x=x+1;
K_TX_nR = params(x); x=x+1;
K_Q = params(x); x=x+1;
h_Q = params(x); x=x+1;
m_deg = params(x); x=x+1;
kb_TL = params(x); x=x+1;
ku_TL = params(x); x=x+1;
v_TL = params(x); x=x+1;
K_TL = params(x); x=x+1;
prom_vec = params(x:x+n-1)'; x=x+n;
z_vec = params(x:x+n-2);


%% VARIABLES

x = 1;

subpop = y(x:x+n-1); x=x+n;
e = y(x:x+n-1); x=x+n;
m_R = y(x:x+n-1); x=x+n;
m_C = y(x:x+n-1); x=x+n;
m_Q = y(x:x+n-1); x=x+n;
m_H = y(x:x+n-1); x=x+n;
TL_R = y(x:x+n-1); x=x+n;
TL_C = y(x:x+n-1); x=x+n;
TL_Q = y(x:x+n-1); x=x+n;
TL_H = y(x:x+n-1); x=x+n;
R = y(x:x+n-1); x=x+n;
C = y(x:x+n-1); x=x+n;
Q = y(x:x+n-1); x=x+n;
H = y(x:x+n-1);



%% CALCULATIONS

% Building blocks
epsilon = nq*C*v_e*nut ./ (K_e + nut);

% Transcription
IQ = (1 ./ (1 + (Q/K_Q).^h_Q));

w_R = (v_TX_R * e) ./ (K_TX_R + e);
w_C = (v_TX_C * e) ./ (K_TX_nR + e);
w_Q = (v_TX_Q * e) ./ (K_TX_nR + e) .* IQ;
w_H = (prom_vec .* e) ./ (K_TX_nR + e);

% Translation
TL_rate = (v_TL * e) ./ (K_TL + e);
gamma_R = TL_R .* TL_rate / n_R;
gamma_C = TL_C .* TL_rate / n_C;
gamma_Q = TL_Q .* TL_rate / n_Q;
gamma_H = TL_H .* TL_rate ./ n_H;

% Summations
m_all = m_R + m_C + m_Q + m_H;
TL_all = TL_R + TL_C + TL_Q + TL_H;
gamma_all = gamma_R + gamma_C + gamma_Q + gamma_H;

% Growth and Dilution
GR = TL_rate .* TL_all / mass;
buffer = sum(subpop) - N;



%% ODEs

dydt(length(y),1) = 0; % Setup ODE structure
x = 1;

for i = 1:n
    
    if i == 1 % First subpop
    dydt(i) = subpop(i)*GR(i)*(1-sum(z_vec))...
            - subpop(i)*buffer;
    
    elseif i < n % Intermediate cells
    dydt(i) = sum(subpop(1:i-1).*GR(1:i-1)*z_vec(i-1))... % Can't use for subpop(1) due to i-1
            + subpop(i)*GR(i)*(1-sum(z_vec(i:end)))... % Can't use for subpop(end) due to z_vec
            - subpop(i)*buffer;
        
    else % Last subpop
    dydt(i) = sum(subpop(1:i-1).*GR(1:i-1)*z_vec(i-1))...
            + subpop(i)*GR(i)...
            - subpop(i)*buffer;
    x=x+n;
    end
end

% e
dydt(x:x+n-1) = epsilon - (gamma_R*n_R + gamma_C*n_C + gamma_Q*n_Q + gamma_H.*n_H)...
                - GR.*e; x=x+n;
% m_R
dydt(x:x+n-1) = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR+m_deg).*m_R; x=x+n;
% m_C
dydt(x:x+n-1) = w_C - kb_TL*m_C.*R + ku_TL*TL_C + gamma_C - (GR+m_deg).*m_C; x=x+n;
% m_Q
dydt(x:x+n-1) = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR+m_deg).*m_Q; x=x+n;
% m_H
dydt(x:x+n-1) = w_H - kb_TL.*m_H.*R + ku_TL.*TL_H + gamma_H - (GR+m_deg).*m_H; x=x+n;
% TL_R
dydt(x:x+n-1) = kb_TL*m_R.*R - ku_TL*TL_R - gamma_R - GR.*TL_R; x=x+n;
% TL_C
dydt(x:x+n-1) = kb_TL*m_C.*R - ku_TL*TL_C - gamma_C - GR.*TL_C; x=x+n;
% TL_Q
dydt(x:x+n-1) = kb_TL*m_Q.*R - ku_TL*TL_Q - gamma_Q - GR.*TL_Q; x=x+n;
% TL_H
dydt(x:x+n-1) = kb_TL.*m_H.*R - ku_TL.*TL_H - gamma_H - GR.*TL_H; x=x+n;
% R
dydt(x:x+n-1) = gamma_R + ku_TL*(TL_all-TL_H) - kb_TL*(m_all-m_H).*R + gamma_all...
                + ku_TL.*TL_H - kb_TL.*m_H.*R...
                - GR.*R; x=x+n;
% C
dydt(x:x+n-1) = gamma_C - GR.*C; x=x+n;
% Q
dydt(x:x+n-1) = gamma_Q - GR.*Q; x=x+n;
% H
dydt(x:x+n-1) = gamma_H - GR.*H;