%%% Additional script that simulates the ODEs for a one-gene system. For
%%% more contextual info, see "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

function [dydt,y] = f_ODEs_1gene(t, y, params, subpops_up, z_values_up, z_values_dn)

%% Unpack components

%%% Parameters
x = 1; % Initialise for automatic numbering
n = params(x); x=x+1;
N = params(x); x=x+1;
mass = params(x); x=x+1;
nut = params(x); x=x+1;
nq = params(x); x=x+1;
v_e = params(x); x=x+1;
K_e = params(x); x=x+1;
n_R = params(x); x=x+1;
n_Z = params(x); x=x+1;
n_Q = params(x); x=x+1;
n_H = params(x); x=x+1;
v_TX_R = params(x); x=x+1;
v_TX_Z = params(x); x=x+1;
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
p_deg = params(x); x=x+1;
alpha_vec = params(x:x+n-1)'; x=x+length(alpha_vec); % Must be column vector
beta_vec = params(x:x+n-1)';

%%% Variables
x = 1; % Initialise for automatic numbering
subpop = y(x:x+n-1); x=x+n;
e = y(x:x+n-1); x=x+n;
m_R = y(x:x+n-1); x=x+n;
m_Z = y(x:x+n-1); x=x+n;
m_Q = y(x:x+n-1); x=x+n;
m_H = y(x:x+n-1); x=x+n;
TL_R = y(x:x+n-1); x=x+n;
TL_Z = y(x:x+n-1); x=x+n;
TL_Q = y(x:x+n-1); x=x+n;
TL_H = y(x:x+n-1); x=x+n;
R = y(x:x+n-1); x=x+n;
Z = y(x:x+n-1); x=x+n;
Q = y(x:x+n-1); x=x+n;
H = y(x:x+n-1);

%% Calculate rates

%%% Energy
epsilon = nq*Z*v_e*nut ./ (K_e + nut);

%%% Autoregulation of m_Q
IQ = (1 ./ (1 + (Q/K_Q).^h_Q)); 

%%% Transcription regulation for synthetic genes. Default is no repression,
%%% but can modify this akin to "Q regultion" above to fit your choen
%%% regulatory network topology.
IH = 1; % No repression

%%% Transcription rates
w_R = (v_TX_R * e) ./ (K_TX_R + e);
w_Z = (v_TX_Z * e) ./ (K_TX_nR + e);
w_Q = (v_TX_Q * e) ./ (K_TX_nR + e) .* IQ;
w_H = (alpha_vec .* e) ./ (K_TX_nR + e) .* IH;

%%% Translation
TL_rate = (v_TL * e) ./ (K_TL + e);
gamma_R = TL_R .* TL_rate / n_R;
gamma_Z = TL_Z .* TL_rate / n_Z;
gamma_Q = TL_Q .* TL_rate / n_Q;
gamma_H = TL_H .* TL_rate ./ n_H;

%%% Useful summations
m_all = m_R + m_Z + m_Q + m_H;
TL_all = TL_R + TL_Z + TL_Q + TL_H;
gamma_all = gamma_R + gamma_Z + gamma_Q + gamma_H;

%%% Growth rate
GR = TL_rate .* TL_all / mass;

%%% Dilution. To ensure that the numerical implementation follows the
%%% mathematical formulation from Equation (1), it can be written without
%%% the if/else structure. If using the if/else structure (i.e. the
%%% mathematically correct formulation), then abs/rel tolerances usually
%%% need to be lowered for the simulations to work.
buffer = sum(subpop) - N;
% if sum(subpop) > N
%     buffer = sum(subpop) - N;
% else
%     buffer = 0;
% end


%% ODEs

%%% For explanations of the 'number of cells in each state', see
%%% Section S1.2. For the cell variables, see Section S1.1.

dydt(length(y),1) = 0; % Initialise to store values
x = 1; % Initialise for automatic numbering

%%% Number of cells in each state
for i = 1:n
    dydt(i) = sum(subpop(subpops_up{i})'.*GR(subpops_up{i})'.*z_values_up{i})... % Division from upstream states
            + subpop(i)*GR(i)*(1-sum(z_values_dn{i}))... % Normal 'self-division'
            - subpop(i)*buffer; % Dilution
end
x=x+n;

%%% e
dydt(x:x+n-1) = epsilon - (gamma_R*n_R + gamma_Z*n_Z + gamma_Q*n_Q...
              + gamma_H.*n_H)...
              - GR.*e; x=x+n;
%%% m_R
dydt(x:x+n-1) = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR+m_deg).*m_R; x=x+n;
%%% m_Z
dydt(x:x+n-1) = w_Z - kb_TL*m_Z.*R + ku_TL*TL_Z + gamma_Z - (GR+m_deg).*m_Z; x=x+n;
%%% m_Q
dydt(x:x+n-1) = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR+m_deg).*m_Q; x=x+n;
%%% m_H
dydt(x:x+n-1) = w_H - beta_vec.*m_H.*R + ku_TL*TL_H + gamma_H - (GR+m_deg).*m_H; x=x+n;
%%% TL_R
dydt(x:x+n-1) = kb_TL*m_R.*R - ku_TL*TL_R - gamma_R - GR.*TL_R; x=x+n;
%%% TL_Z
dydt(x:x+n-1) = kb_TL*m_Z.*R - ku_TL*TL_Z - gamma_Z - GR.*TL_Z; x=x+n;
%%% TL_Q
dydt(x:x+n-1) = kb_TL*m_Q.*R - ku_TL*TL_Q - gamma_Q - GR.*TL_Q; x=x+n;
%%% TL_H
dydt(x:x+n-1) = beta_vec.*m_H.*R - ku_TL*TL_H - gamma_H - GR.*TL_H; x=x+n;
%%% R
dydt(x:x+n-1) = gamma_R + ku_TL*(TL_all-TL_H) - kb_TL*(m_all-m_H).*R + gamma_all...
              + ku_TL*TL_H - beta_vec.*m_H.*R...
              - GR.*R; x=x+n;
%%% Z
dydt(x:x+n-1) = gamma_Z - GR.*Z; x=x+n;
%%% Q
dydt(x:x+n-1) = gamma_Q - GR.*Q; x=x+n;
%%% H
dydt(x:x+n-1) = gamma_H - (GR+p_deg).*H;