%% Care

% Any parameter vectors I use in ODEs must be column vectors, not row
% vectors, e.g. prom_vec and RBS_vec. E.g. prom_vec*plasmid*P =
% [a,b,c]*1*[x;y;z] = single value (3*3 matrix multplication).

function [dydt,y] = MWCM2_odes_repressilator(t, y, params, subpops_up, z_values_up, z_values_dn)

%% Unpack params, populate variables

%%%%%%%%%% PARAMETERS

x = 1;

N = params(x); x=x+1;
n = params(x); x=x+1;
mass = params(x); x=x+1;
nut = params(x); x=x+1;
nq = params(x); x=x+1;
v_e = params(x); x=x+1;
K_e = params(x); x=x+1;
n_R = params(x); x=x+1;
n_C = params(x); x=x+1;
n_Q = params(x); x=x+1;
n_HA = params(x); x=x+1;
n_HB = params(x); x=x+1;
n_HC = params(x); x=x+1;
v_TX_R = params(x); x=x+1;
v_TX_C = params(x); x=x+1;
v_TX_Q = params(x); x=x+1;
K_TX_R = params(x); x=x+1;
K_TX_nR = params(x); x=x+1;
K_Q = params(x); x=x+1;
h_Q = params(x); x=x+1;
K_H = params(x); x=x+1;
h_H = params(x); x=x+1;
m_deg = params(x); x=x+1;
kb_TL = params(x); x=x+1;
ku_TL = params(x); x=x+1;
v_TL = params(x); x=x+1;
K_TL = params(x); x=x+1;
p_deg = params(x); x=x+1;

prom_plusA_vec = params(x:x+n-1)'; x=x+length(prom_plusA_vec); % RowVec -> ColVec etc.
RBS_plusA_vec = params(x:x+n-1)'; x=x+length(RBS_plusA_vec);
RBS_minusA_vec = params(x:x+n-1)'; x=x+length(RBS_plusA_vec);
CDSA_vec = params(x:x+n-1)'; x=x+length(CDSA_vec);
prom_plusB_vec = params(x:x+n-1)'; x=x+length(prom_plusB_vec); % RowVec -> ColVec etc.
RBS_plusB_vec = params(x:x+n-1)'; x=x+length(RBS_plusB_vec);
RBS_minusB_vec = params(x:x+n-1)'; x=x+length(RBS_plusB_vec);
CDSB_vec = params(x:x+n-1)'; x=x+length(CDSB_vec);
prom_plusC_vec = params(x:x+n-1)'; x=x+length(prom_plusC_vec); % RowVec -> ColVec etc.
RBS_plusC_vec = params(x:x+n-1)'; x=x+length(RBS_plusC_vec);
RBS_minusC_vec = params(x:x+n-1)'; x=x+length(RBS_plusC_vec);
CDSC_vec = params(x:x+n-1)'; % x=x+length(CDSC_vec);



%% VARIABLES

x = 1;

subpop = y(x:x+n-1); x=x+n;
e = y(x:x+n-1); x=x+n;
m_R = y(x:x+n-1); x=x+n;
m_C = y(x:x+n-1); x=x+n;
m_Q = y(x:x+n-1); x=x+n;
m_HA = y(x:x+n-1); x=x+n;
m_HB = y(x:x+n-1); x=x+n;
m_HC = y(x:x+n-1); x=x+n;
TL_R = y(x:x+n-1); x=x+n;
TL_C = y(x:x+n-1); x=x+n;
TL_Q = y(x:x+n-1); x=x+n;
TL_HA = y(x:x+n-1); x=x+n;
TL_HB = y(x:x+n-1); x=x+n;
TL_HC = y(x:x+n-1); x=x+n;
R = y(x:x+n-1); x=x+n;
C = y(x:x+n-1); x=x+n;
Q = y(x:x+n-1); x=x+n;
HA = y(x:x+n-1); x=x+n;
HB = y(x:x+n-1); x=x+n;
HC = y(x:x+n-1);



%% CALCULATIONS

% Building blocks
epsilon = nq*C*v_e*nut ./ (K_e + nut);

% Transcription
IQ = (1 ./ (1 + (Q/K_Q).^h_Q));         % AutoInhibition of m_Q
IHA = (1 ./ (1 + (HC/K_H).^h_H));       % AutoInhibition of m_HA via HC
IHB = (1 ./ (1 + (HA/K_H).^h_H));       % AutoInhibition of m_HB via HA
IHC = (1 ./ (1 + (HB/K_H).^h_H));       % AutoInhibition of m_HC via HB

w_R = (v_TX_R * e) ./ (K_TX_R + e);
w_C = (v_TX_C * e) ./ (K_TX_nR + e);
w_Q = (v_TX_Q * e) ./ (K_TX_nR + e) .* IQ;
w_HA = (prom_plusA_vec .* e) ./ (K_TX_nR + e) .* IHA;
w_HB = (prom_plusB_vec .* e) ./ (K_TX_nR + e) .* IHB;
w_HC = (prom_plusC_vec .* e) ./ (K_TX_nR + e) .* IHC;

% Translation
TL_rate = (v_TL * e) ./ (K_TL + e);
gamma_R = TL_R .* TL_rate / n_R;
gamma_C = TL_C .* TL_rate / n_C;
gamma_Q = TL_Q .* TL_rate / n_Q;
gamma_HA = TL_HA .* TL_rate ./ CDSA_vec;
gamma_HB = TL_HB .* TL_rate ./ CDSB_vec;
gamma_HC = TL_HC .* TL_rate ./ CDSC_vec;

% Summations
w_all = w_R + w_C + w_Q + w_HA + w_HB + w_HC;
m_all = m_R + m_C + m_Q + m_HA + m_HB + m_HC;
TL_all = TL_R + TL_C + TL_Q + TL_HA + TL_HB + TL_HC;
gamma_all = gamma_R + gamma_C + gamma_Q + gamma_HA + gamma_HB + gamma_HC;

% Growth and Dilution
GR = TL_rate .* TL_all / mass;
buffer = sum(subpop) - N;



%% ODEs

dydt(length(y),1) = 0; % Setup ODE structure
x = 1;

for i = 1:n    
    
    % One eqn form works for all subpops, as when no cells exist, values
    % are just not included in the calculation
    dydt(i) = sum(subpop(subpops_up{i})'.*GR(subpops_up{i})'.*z_values_up{i})...
            + subpop(i)*GR(i)*(1-sum(z_values_dn{i}))...
            - subpop(i)*buffer;
end

x=x+n;

% e
dydt(x:x+n-1) = epsilon - (gamma_R*n_R + gamma_C*n_C + gamma_Q*n_Q...
              + gamma_HA.*CDSA_vec + gamma_HB.*CDSB_vec + gamma_HC.*CDSC_vec)...
              - GR.*e; x=x+n;
% m_R
dydt(x:x+n-1) = w_R - kb_TL*m_R.*R + ku_TL*TL_R + gamma_R - (GR+m_deg).*m_R; x=x+n;
% m_C
dydt(x:x+n-1) = w_C - kb_TL*m_C.*R + ku_TL*TL_C + gamma_C - (GR+m_deg).*m_C; x=x+n;
% m_Q
dydt(x:x+n-1) = w_Q - kb_TL*m_Q.*R + ku_TL*TL_Q + gamma_Q - (GR+m_deg).*m_Q; x=x+n;
% m_HA
dydt(x:x+n-1) = w_HA - RBS_plusA_vec.*m_HA.*R + RBS_minusA_vec.*TL_HA + gamma_HA - (GR+m_deg).*m_HA; x=x+n;
% m_HB
dydt(x:x+n-1) = w_HB - RBS_plusB_vec.*m_HB.*R + RBS_minusB_vec.*TL_HB + gamma_HB - (GR+m_deg).*m_HB; x=x+n;
% m_HC
dydt(x:x+n-1) = w_HC - RBS_plusC_vec.*m_HC.*R + RBS_minusC_vec.*TL_HC + gamma_HC - (GR+m_deg).*m_HC; x=x+n;
% TL_R
dydt(x:x+n-1) = kb_TL*m_R.*R - ku_TL*TL_R - gamma_R - GR.*TL_R; x=x+n;
% TL_C
dydt(x:x+n-1) = kb_TL*m_C.*R - ku_TL*TL_C - gamma_C - GR.*TL_C; x=x+n;
% TL_Q
dydt(x:x+n-1) = kb_TL*m_Q.*R - ku_TL*TL_Q - gamma_Q - GR.*TL_Q; x=x+n;
% TL_HA
dydt(x:x+n-1) = RBS_plusA_vec.*m_HA.*R - RBS_minusA_vec.*TL_HA - gamma_HA - GR.*TL_HA; x=x+n;
% TL_HB
dydt(x:x+n-1) = RBS_plusB_vec.*m_HB.*R - RBS_minusB_vec.*TL_HB - gamma_HB - GR.*TL_HB; x=x+n;
% TL_HC
dydt(x:x+n-1) = RBS_plusC_vec.*m_HC.*R - RBS_minusC_vec.*TL_HC - gamma_HC - GR.*TL_HC; x=x+n;
% R
dydt(x:x+n-1) = gamma_R + ku_TL*(TL_all-(TL_HA+TL_HB+TL_HC)) - kb_TL*(m_all-(m_HA+m_HB+m_HC)).*R + gamma_all...
              + RBS_minusA_vec.*TL_HA - RBS_plusA_vec.*m_HA.*R...
              + RBS_minusB_vec.*TL_HB - RBS_plusB_vec.*m_HB.*R...
              + RBS_minusC_vec.*TL_HC - RBS_plusC_vec.*m_HC.*R...
              - GR.*R; x=x+n;
% C
dydt(x:x+n-1) = gamma_C - GR.*C; x=x+n;
% Q
dydt(x:x+n-1) = gamma_Q - GR.*Q; x=x+n;

%%% Protein degradation
% HA
dydt(x:x+n-1) = gamma_HA - (GR+p_deg).*HA; x=x+n;
% HB
dydt(x:x+n-1) = gamma_HB - (GR+p_deg).*HB; x=x+n;
% HC
dydt(x:x+n-1) = gamma_HC - (GR+p_deg).*HC;

%%% No protein degradtation
% % HA
% dydt(x:x+n-1) = gamma_HA - (GR).*HA; x=x+n;
% % HB
% dydt(x:x+n-1) = gamma_HB - (GR).*HB; x=x+n;
% % HC
% dydt(x:x+n-1) = gamma_HC - (GR).*HC;