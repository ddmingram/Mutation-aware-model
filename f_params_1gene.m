%%% Additional script that initialises parameters for the mutation model
%%% when expressing one synthetic gene. For more contextual info, see
%%% "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

%%% Optional but useful
close; clear; clc; % Clear everything
tic % Starts timer

%%% Normal parameters - these will be present for all standard simulations
N = 1e9;                % Cell capacity of turbidostat
mass = 10^8;            % Cell mass (aa)
nut = 1e4;              % Internal nutrient quantity (molecs)
nq = 1;                 % Nutrient quality
v_e = 38700;            % Max rate of energy metabolism (/h)
K_e = 1000;             % Half-saturation constant for energy metabolism (molecs)
n_R = 7459;             % Length of R-proteins (aa/molecs)
n_Z = 300;              % Length of Z-proteins (aa/molecs)
n_Q = 300;              % Length of Q-proteins (aa/molecs)
n_H = 300;              % Length of H-proteins (aa/molecs)
v_TX_R = 55800;         % Max rate of R-transcription (molecs/h/cell)
v_TX_Z = 248.4;         % Max rate of Z-transcription (molecs/h/cell)
v_TX_Q = 56940;         % Max rate of Q-transcription (molecs/h/cell)
K_TX_R = 427;           % Half-saturation constant for R-transcription (molecs/cell)
K_TX_nR = 4.38;         % Half-saturation constant for non-R-transcription (molecs/cell)
K_Q = 152000;           % Half-saturation constant for Q-inhibition (molecs/cell)
h_Q = 4;                % Hill coefficient for Q-inhibition
m_deg = 60*(log(2)/7);  % mRNA degradation rate (/h) (= deg. time of 7min)
kb_TL = 60;             % mRNA-ribosome binding rate for R/Z/Q (cell/h/molecs)
ku_TL = 60;             % mRNA-ribosome unbinding rate for R/Z/Q (/h)
v_TL = 72000;           % Max rate of translation (aa/h/molecs)
K_TL = 7;               % Half-saturation constant for translation (molecs/cell)
p_deg = 0;              % Protein degradation rate (/h)

%%% Consolidate into one variable
base_params = [N, mass, nut, nq, v_e, K_e, n_R, n_Z, n_Q, n_H, v_TX_R,...
               v_TX_Z, v_TX_Q, K_TX_R, K_TX_nR, K_Q, h_Q, m_deg, kb_TL,...
               ku_TL, v_TL, K_TL, p_deg];