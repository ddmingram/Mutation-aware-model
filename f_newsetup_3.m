%%% A script that initialises conditions to run the mutation model,
%%% specific for *three gene constructs*. These may be on the same or
%%% separate plasmids.

%% Set-up

close; clear;

tic % Start timer

set(0,'DefaultLineLineWidth',4);
set(0,'defaultAxesFontSize',32);

E_yellow = [224 199 5]/255;
I_orange = [224 135 0]/255;
M_red = [255, 51, 51]/255;
LightGrey = [230 230 230]/255;
MediumGrey = [175 175 175]/255;
DarkGrey = [51 51 51]/255;
MATLABblue = [3 115 189]/255;
MATLABblueFaded = [93 169 221]/255;
MATLABorange = [217 87 30]/255;
MATLABorangeFaded = [234 148 111]/255;
MATLAByellow = [237 177 32]/255;
MATLABpurple = [126 47 142]/255;



%% Base parameters

a = 1;                  % Change to 60 to convert to /min

N = 1e9;
n = 2;                  % Value will be updated in core script
mass = 10^8;
nut = 1e4;              % Extracellular nutrient
nq = 1;                 % Nutrient Quality
v_e = 38700/a;          % /h ...= 1/((1/trans)+(1/met))
K_e = 500;             % molecs/cell
n_R = 7459;
n_C = 300;
n_Q = 300;
n_HA = 300;
n_HB = 300;
n_HC = 300;
v_TX_R = 55800/a;       % molecs/h/cell
v_TX_C = 248.4/a;         % molecs/h/cell
v_TX_Q = 56940/a;       % molecs/h/cell
K_TX_R = 427;           % molecs/cell
K_TX_nR = 4.38;            % molecs/cell
K_Q = 152000;           % molecs/cell
h_Q = 4;
K_H = 100;
h_H = 2;
m_deg = 60*(log(2)/2)/a;       % /h
kb_TL = 60/a;           % cell/h/molecs
ku_TL = 60/a;           % /h
v_TL = 72000/a;         % aa/h (=20 aa/s)
K_TL = 7;               % molecs/cell (=v_TL/K_P = 72000/(180*60)) 
p_deg = 0.4*60*(log(2)/4);              % Protein deg' 0 by default
% p_deg = 0;

base_params = [N, n, mass,...
               nut, nq, v_e, K_e, n_R, n_C, n_Q, n_HA, n_HB, n_HC,...
               v_TX_R, v_TX_C, v_TX_Q, K_TX_R, K_TX_nR, K_Q, h_Q, K_H, h_H, m_deg,...
               kb_TL, ku_TL, v_TL, K_TL, p_deg];


%% ODE options

tspan = [0, 120];

Opt_1 = odeset('RelTol',1e-6,'AbsTol',1e-9);
% Opt_1 = odeset('RelTol',1e-15,'AbsTol',1e-18);
% Default: RelTol: 1e-3, AbsTol: 1e-6
% Increase from default to avoid weird plitting

Opt_2 = odeset('Events', @f_H_terminate);
% Terminate ODE early if H falls below threshold!

Opt_all = odeset(Opt_1, Opt_2);