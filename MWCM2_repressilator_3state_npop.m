% Run the set-up file
f_newsetup_3;

ABC_A = [95 126 196]/255;
ABC_B = [252 141 98]/255;
ABC_C = [88 178 150]/255;



%% 1. Define parts and probabilities to consider

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parts_vals_mat = [d1s1  d2s1  ...
%                   d1s2  d2s2  ...
%                     :     :      ]
%
% z_mat = [d1,->s1  d2,->s1  ...    -> = going to state x
%          d1,->s2  d2,->s2  ...  ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Record all values that I might use
% parts_mat_all = [[69; 6.9; 0],...       % prom_A
%                  [69; 1; 0],...         % RBS_A
%                  [1; 0.2; 0.1],...      % CDS_A
%                  [69; 6.9; 0],...       % prom_B
%                  [69; 1; 0],...         % RBS_B
%                  [1; 0.2; 0.1],...      % CDS_B
%                  [69; 6.9; 0],...       % prom_C
%                  [69; 1; 0],...         % RBS_C
%                  [1; 0.2; 0.1]];        % CDS_C

%%% For repressilator comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To plot in /min, set induction range to 60eX-->60eY. Then divide values
% by 60 when plotting, and set x-range to 1eX-->1eY

trials = 30;
Amplitude_store = zeros(1,trials); 
Period_store = zeros(1,trials);
induc_vec = logspace(0,4,trials);
% induc_vec = logspace(0,4,trials)*60;

% for k = 1:trials
for k = 1
    
induc_vec = 60*69;
    
% induc_vec = 6000;
% parts_mat_all = [[induc_vec(k); 60*34.5; 0],...       % prom_A
%                  [60; 0.5; 0],...         % RBS_A
%                  [1; 0.2; 0.1],...      % CDS_A
%                  [induc_vec(k); 60*34.5; 0],...       % prom_B
%                  [60; 0.5; 0],...         % RBS_B
%                  [1; 0.2; 0.1],...      % CDS_B
%                  [induc_vec(k); 60*34.5; 0],...       % prom_C
%                  [60; 0.5; 0],...         % RBS_C
%                  [1; 0.2; 0.1]];        % CDS_C

part_WT = 1e5;
part_int = 1e4;

parts_mat_all = [[part_WT; part_int; 0],...       % prom_A
                 [60; 6; 0],...         % RBS_A
                 [1; 0.2; 0.1],...      % CDS_A
                 [part_WT; part_int; 0],...       % prom_B
                 [60; 6; 0],...         % RBS_B
                 [1; 0.2; 0.1],...      % CDS_B
                 [part_WT; part_int; 0],...       % prom_C
                 [60; 6; 0],...         % RBS_C
                 [1; 0.2; 0.1]];        % CDS_C
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
z_severe = 1e-3;
z_partial = 5e-2;

z_mat_all = [[z_severe; z_partial],...           % prom_A
             [z_severe; z_partial],...           % Etc.
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial],...
             [z_severe; z_partial]];

% z_mat_all = [[1e-3; 1e-1],...           % prom_A
%              [1e-3; 1e-1],...           % Etc.
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1],...
%              [1e-3; 1e-1]];

%%% Essential params
dim_vec = [1,0,0, 1,0,0, 1,0,0];  % Which dims are active? idx=1 -> active.
d = sum(dim_vec);   % How many dims are there?
s = 3;              % How many states do we want?
n = s^d;



%% 2. Run functions to get all info needed...

%%% Get structure of the subpop network
MyPop = f_grid_structure(s,d);

%%% Order and distribute part values to params for each subpop
[z_mat, parts_mat,...
 prom_plusA_vec,...
 prom_minusA_vec,...
 RBS_plusA_vec,...
 RBS_minusA_vec,...
 CDSA_vec,...
 prom_plusB_vec,...
 prom_minusB_vec,...
 RBS_plusB_vec,...
 RBS_minusB_vec,...
 CDSB_vec,...
 prom_plusC_vec,...
 prom_minusC_vec,...
 RBS_plusC_vec,...
 RBS_minusC_vec,...
 CDSC_vec] = f_PartParameters_3(dim_vec, s, d, n,...
                                parts_mat_all,...
                                z_mat_all,...
                                n_HA, n_HB, n_HC);

%%% Order and distribute part values to params for each subpop
[subpops_up,...
 z_values_up,...
 subpops_dn,...
 z_values_dn] = f_popODEinfo(z_mat, MyPop, s, d, n);


                            
%% 3. Collate parameters

%%% Collect parameters, replacing n from base_params
base_params(2) = n;

%% Repressilator vary info
    
params = [base_params,...
          prom_plusA_vec,...
          RBS_plusA_vec, RBS_minusA_vec,...
          CDSA_vec,...
          prom_plusB_vec,...
          RBS_plusB_vec, RBS_minusB_vec,...
          CDSB_vec,...
          prom_plusC_vec,...
          RBS_plusC_vec, RBS_minusC_vec,...
          CDSC_vec];

      

%% 4. Initialise and store variables

subpop_0 = [N, zeros(1,n-1)]; % Start with all first population (E)

%%% Assign SS values - if H starts at 0
e_0 = [9.94582659197826,10.7810881240108,9.69973920909995,10.7047883481363,12.0024842493255,10.8941407868587,9.69973920910001,9.69974321083094,9.69974365539743,9.76051402483775,10.8883388269577,9.69974321083041,10.2590291805889,12.0452720578585,10.8946469149369,10.8941407868587,10.8946469149370,10.8947031362534,9.69973920909992,10.8941407868587,9.69974365539720,9.69974321083045,10.8946469149369,10.8947031362534,9.69974365539706,10.8947031362535,23.0502435544397]; 
m_R_0 = [65.2697549460658,47.2232971701486,58.8736101791348,65.1075415644158,39.7158129127373,51.1032910617358,58.8736101791381,58.8735634924139,58.8735583059472,66.8900995072828,49.3356438797402,58.8735634924148,82.0196434478997,51.1677364468182,51.1005050043859,51.1032910617360,51.1005050043862,51.1001955404294,58.8736101791420,51.1032910617357,58.8735583059506,58.8735634924225,51.1005050043858,51.1001955404294,58.8735583059573,51.1001955404295,12.1015042756136];
m_C_0 = [8.43479217465313,7.31810075391791,8.14377167604025,8.60840610868366,6.55276431576519,7.55325313543652,8.14377167604062,8.14376600307832,8.14376537287442,8.23059984989335,7.36624607793049,8.14376600307863,8.69016901815503,7.60359311708629,7.55304828588574,7.55325313543652,7.55304828588574,7.55302553127261,8.14377167604112,7.55325313543652,8.14376537287492,8.14376600307953,7.55304828588572,7.55302553127262,8.14376537287578,7.55302553127263,2.50884469718633];

m_Q_0 = [1067.91272565901,653.634748836675,1544.29482410476,670.707563760963,410.994093690680,565.987539244195,1544.29482410483,1544.25118058690,1544.24633236513,805.258262645734,595.168265519489,1544.25118058696,679.962010595144,464.271648262520,565.910046975431,565.987539244196,565.910046975432,565.901441213382,1544.29482410499,565.987539244195,1544.24633236522,1544.25118058719,565.910046975430,565.901441213381,1544.24633236545,565.901441213382,96.5367857954868];
m_HA_0 = [6.30414390029838,68.7425616979676,-1.50000000000000e-323,537.171447243256,58.8743246546454,-4.90000000000000e-324,0.328138920375014,0.0328052202558561,4.00000000000000e-323,2774.00331158762,242.356232514732,1.50000000000000e-323,2228.14743817317,4.12098989125331,-1.50000000000000e-323,0.101209432208117,0.0101151144802782,1.50000000000000e-323,3278.49101289900,304.076213181824,1.50000000000000e-323,3278.48872909804,304.067966420520,4.90000000000000e-324,3278.48847539282,304.067050373294,5.62959709939780e-162];
m_HB_0 = [1.31276681442576,0.566930584772558,3278.49101289865,4.97751729217057,0.811073462173470,304.076213181824,-1.50000000000000e-323,4.00000000000000e-323,0,0.524225054378654,0.279346057729880,3278.48872909768,0.772582624703915,0.300609625690701,304.067966420521,1.00000000000000e-323,4.90000000000000e-324,4.90000000000000e-324,0.328138920379140,0.101209432208117,3278.48847539248,0.0328052202568932,0.0101151144802782,304.067050373294,-1.50000000000000e-323,-4.90000000000000e-324,3.31902464938638e-162];

m_HC_0 = [1542.77709148554,40.1101020923173,0.328138920372218,4.19458662964150,19.0250871889118,0.101209432208117,3278.49101289880,3278.48872909755,3278.48847539228,1.01574286519183,5.21910690867773,0.0328052202562030,2.81453446069162,248.368993833309,0.0101151144802782,304.076213181825,304.067966420521,304.067050373293,-1.50000000000000e-323,1.50000000000000e-323,-4.00000000000000e-323,-1.50000000000000e-323,1.50000000000000e-323,4.90000000000000e-324,-1.50000000000000e-323,-2.00000000000000e-323,8.44465157489357e-165];
TL_R_0 = [69.7119319196450,233.955315615268,17.0219473815796,230.332532272373,431.756752869371,312.405059010307,17.0219473815085,17.0242388352992,17.0244933976484,86.5618131042059,274.722635047429,17.0242388352101,176.456082858604,457.186403298416,312.488000140593,312.405059010307,312.488000140594,312.497212365686,17.0219473814031,312.405059010306,17.0244933975342,17.0242388350342,312.488000140593,312.497212365686,17.0244933973582,312.497212365686,1187.57007996983];
TL_C_0 = [2.89536785797530,11.7961103721956,0.776287951361163,9.77108404641987,22.5884402389200,14.8938854028331,0.776287951357907,0.776392552436314,0.776404172788805,3.43247508990180,13.2687509043744,0.776392552432282,5.91497058372594,21.4135636759941,14.8981009090429,14.8938854028330,14.8981009090428,14.8985691337136,0.776287951353102,14.8938854028331,0.776404172783606,0.776392552424240,14.8981009090429,14.8985691337136,0.776404172775574,14.8985691337136,69.5589263082164];

TL_Q_0 = [366.661094337474,1056.10187424911,147.206664551877,761.763627595665,1418.20854438615,1116.04276962278,147.206664551259,147.222441711299,147.224194382568,334.368698879147,1073.05333024759,147.222441710534,459.459909707485,1306.63319907973,1116.23607663624,1116.04276962278,1116.23607663624,1116.25754604926,147.206664550354,1116.04276962278,147.224194381583,147.222441709015,1116.23607663624,1116.25754604925,147.224194380065,1116.25754604925,2676.52883285729];
TL_HA_0 = [2.32469562746880,115.273822385255,0,596.307094303088,208.359613916458,-1.00000000000000e-323,0.0312791542287640,0.00312751234233050,4.90000000000000e-324,1149.06332000600,442.285956050501,0,1500.00068255417,11.9507227368372,-3.00000000000000e-323,0.199569861881292,0.0199516791803527,3.00000000000000e-323,312.515278322505,599.592810097950,0,312.557388254525,599.762516467106,1.00000000000000e-323,312.562066334772,599.781366091529,1.57935542590910e-160];
TL_HB_0 = [0.431117772410498,0.910016981529819,312.515278325750,5.91717226348334,2.87544209959963,599.592810097950,0,4.90000000000000e-324,0,0.234096269217419,0.511995166609969,312.557388257763,0.557734945985464,0.814973413941771,599.762516467105,2.00000000000000e-323,1.00000000000000e-323,1.00000000000000e-323,0.0312791542289621,0.199569861881291,312.562066338006,0.00312751234238038,0.0199516791803528,599.781366091530,0,-1.00000000000000e-323,9.31135833734029e-161];

TL_HC_0 = [517.965630049057,62.0184541250063,0.0312791542286298,4.56808670324216,63.0319975943720,0.199569861881291,312.515278324439,312.557388259386,312.562066340099,0.403603588442398,8.99864977467875,0.00312751234234725,1.82113329672192,696.045515308687,0.0199516791803528,599.592810097949,599.762516467103,599.781366091530,0,3.00000000000000e-323,-4.90000000000000e-324,0,3.00000000000000e-323,1.00000000000000e-323,0,-4.00000000000000e-323,2.36910493757123e-163];
R_0 = [1.13791291994864,5.56232789581102,0.317106469876517,3.86361029177353,12.2531897847225,6.80321636988358,0.317106469875173,0.317149479087636,0.317154257074951,1.38276841374408,6.23671160878227,0.317149479085972,2.28371311143081,9.97614574327890,6.80541822011707,6.80321636988355,6.80541822011705,6.80566279646055,0.317106469873190,6.80321636988358,0.317154257072805,0.317149479082653,6.80541822011708,6.80566279646055,0.317154257069489,6.80566279646053,113.797442173045];
C_0 = [1100.66618257758,1753.04546007823,541.852917553901,1900.03832703879,2648.91150499983,2429.90828756700,541.852917551631,541.889469220578,541.893529670707,1790.68258846004,2155.57076915263,541.889469217751,2489.82132503787,3081.37450275180,2430.27422222836,2429.90828756700,2430.27422222836,2430.31486481764,541.852917548282,2429.90828756700,541.893529667072,541.889469212147,2430.27422222836,2430.31486481764,541.893529661469,2430.31486481764,5894.33797183082];

Q_0 = [144429.158293806,173822.346544800,102750.996625015,180282.204711438,196379.965386146,182080.196123289,102750.996625015,102755.095402538,102755.550706705,159512.796090225,178483.315782360,102755.095402538,171161.501986643,194315.799633870,182087.621740015,182080.196123289,182087.621740015,182088.446399159,102750.996625009,182080.196123289,102755.550706705,102755.095402532,182087.621740015,182088.446399159,102755.550706699,182088.446399159,226805.765550571];
HA_0 = [4035.59938552051,6869.73518285832,1.36625349077755e-224,8643.85401368863,8862.06918715539,9.96388029830398e-260,1.00039183552448,0.100026000077101,1.18224839016474e-224,26980.9596371266,13469.7499486544,1.36308323493149e-224,32528.1350802297,2112.78939910973,9.63871910281119e-260,5.76915281002215,0.576756958124309,8.85564717747083e-260,9995.08268746753,17332.9906261259,1.36273128954234e-224,9996.40030785042,17337.7489417124,9.60328812965230e-260,9996.54668470561,17338.2774576021,8.44420942029446e-159];
HB_0 = [84.2507340741262,609.812781123375,9995.08268757133,2091.18951188777,764.560359094757,17332.9906261259,5.47361413615988e-224,4.51425283796816e-224,4.35474867725743e-224,1397.55903834605,517.331395960134,9996.40030795394,810.953013151830,24.2429975037990,17337.7489417124,7.16328281818023e-259,3.76390617740547e-259,3.26104605678458e-259,1.00039183553142,5.76915281002214,9996.54668480909,0.100026000078845,0.576756958124309,17338.2774576021,6.04857363475976e-224,3.95059261870797e-259,4.97842717972410e-159];

HC_0 = [8409.15312473655,608.033496416748,1.00039183551979,184.313952715969,662.862342404688,5.76915281002214,9995.08268752941,9996.40030800607,9996.54668487613,31.1164309785138,89.5824174695456,0.100026000077685,57.4211010232334,14801.3636055060,0.576756958124309,17332.9906261258,17337.7489417124,17338.2774576021,1.89982712231159e-225,1.26115428015800e-260,1.55949250903451e-225,1.89621690788940e-225,1.22780799053711e-260,1.09410221192496e-260,1.89581621345635e-225,1.22416810701164e-260,1.26666980106709e-161];

% ss_store = [e(end,1);...
%          m_R(end,1);...
%          m_C(end,1);...
%          m_Q(end,1);...
%          m_HA(end,1);...
%          m_HB(end,1);...
%          m_HC(end,1);...
%          TL_R(end,1);...
%          TL_C(end,1);...
%          TL_Q(end,1);...
%          TL_HA(end,1);...
%          TL_HB(end,1);...
%          TL_HC(end,1);...
%          R(end,1);...
%          C(end,1);...
%          Q(end,1);...
%          HA(end,1);...
%          HB(end,1);...
%          HC(end,1)];

var = [subpop_0, e_0,...
       m_R_0, m_C_0, m_Q_0, m_HA_0, m_HB_0, m_HC_0,...
       TL_R_0, TL_C_0, TL_Q_0, TL_HA_0, TL_HB_0, TL_HC_0,...
       R_0, C_0, Q_0, HA_0, HB_0, HC_0];

%%% Allocate variables to store results
GR_store = zeros(1,1,n); % Store GR of each subpop!
GR_ratio_store = zeros(1,1,n); % GR ratio between E and M



%% 5. Call ODE solver

[t,y] = ode15s( @(t,y) MWCM2_odes_repressilator(t, y, params, subpops_up, z_values_up, z_values_dn),...
               tspan,...
               var,...
               Opt_1);


           
%% 6. Extract variables

x=1;

subpop = y(:, x:x+n-1); x=x+n;
e = y(:, x:x+n-1); x=x+n;
m_R = y(:, x:x+n-1); x=x+n;
m_C = y(:, x:x+n-1); x=x+n;
m_Q = y(:, x:x+n-1); x=x+n;
m_HA = y(:, x:x+n-1); x=x+n;
m_HB = y(:, x:x+n-1); x=x+n;
m_HC = y(:, x:x+n-1); x=x+n;
TL_R = y(:, x:x+n-1); x=x+n;
TL_C = y(:, x:x+n-1); x=x+n;
TL_Q = y(:, x:x+n-1); x=x+n;
TL_HA = y(:, x:x+n-1); x=x+n;
TL_HB = y(:, x:x+n-1); x=x+n;
TL_HC = y(:, x:x+n-1); x=x+n;
R = y(:, x:x+n-1); x=x+n;
C = y(:, x:x+n-1); x=x+n;
Q = y(:, x:x+n-1); x=x+n;
HA = y(:, x:x+n-1); x=x+n;
HB = y(:, x:x+n-1); x=x+n;
HC = y(:, x:x+n-1);

toc % End timer



%% 7. Calculations

% Building blocks
epsilon = nq*C*v_e*nut ./ (K_e + nut);

% Transcription
IQ = (1 ./ (1 + (Q/K_Q).^h_Q));         % AutoInhibition of m_Q
IHA = (1 ./ (1 + (HC/K_H).^h_H));       % AutoInhibition of m_HA via HC
IHB = (1 ./ (1 + (HA/K_H).^h_H));       % AutoInhibition of m_HB via HA
IHC = (1 ./ (1 + (HB/K_H).^h_H));       % AutoInhibition of m_HC via HB

w_R = (v_TX_R * e) ./ (K_TX_R + e);
w_C = (v_TX_C * e) ./ (K_TX_R + e);
w_Q = (v_TX_Q * e) ./ (K_TX_R + e) .* IQ;
w_HA = (v_TX_R * e) ./ (K_TX_R + e) .* IHA;
w_HB = (v_TX_R * e) ./ (K_TX_R + e) .* IHB;
w_HC = (v_TX_R * e) ./ (K_TX_R + e) .* IHC;

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
buffer = sum(subpop,2) - N;


%% H_avg calcalations

if dim_vec(3) == 1
    % If CDS mutants present, only use non-CDS mutants in calculation
    HA_pop_percell = sum(HA(:,noHmutA_idx).*subpop(:,noHmutA_idx), 2) ./ sum(subpop,2);
    HA_pop_total = sum(HA(:,noHmutA_idx).*subpop(:,noHmutA_idx),2);
else
    HA_pop_percell = sum(HA.*subpop, 2) ./ sum(subpop,2);
    HA_pop_total = sum(HA.*subpop,2);
end

if dim_vec(6) == 1
    % If CDS mutants present, only use non-CDS mutants in calculation
    HB_pop_percell = sum(HB(:,noHmutB_idx).*subpop(:,noHmut_idx), 2) ./ sum(subpop,2);
    HB_pop_total = sum(HB(:,noHmutA_idx).*subpop(:,noHmutB_idx),2);
else
    HB_pop_percell = sum(HB.*subpop, 2) ./ sum(subpop,2);
    HB_pop_total = sum(HB.*subpop,2);
end

if dim_vec(9) == 1
    % If CDS mutants present, only use non-CDS mutants in calculation
    HC_pop_percell = sum(HC(:,noHmutC_idx).*subpop(:,noHmut_idx), 2) ./ sum(subpop,2);
    HC_pop_total = sum(HC(:,noHmutC_idx).*subpop(:,noHmutC_idx),2);
else
    HC_pop_percell = sum(HC.*subpop, 2) ./ sum(subpop,2);
    HC_pop_total = sum(HC.*subpop,2);
end

H_pop_percell = HA_pop_percell + HB_pop_percell + HC_pop_percell;
% H_avg_per_subpop = (HA+HB+HC)/3;

%%% Create a smoothed version of total protein for better calcs
H_pop_percell_S1e3 = movmean(H_pop_percell,1000);

GR_pop_percell = sum(GR.*subpop, 2) ./ sum(subpop,2);
rel_GR_pop_percell = GR_pop_percell ./ GR(:,end);



%% Calculate peaks

% [Amp,Idx] = findpeaks(HA(:,1));
% t_peaks = t(Idx);
% AmpRange = ceil(length(Amp)*0.8); % Only consider last 20%
% tpkRange = ceil(length(t_peaks)*0.8); % Only consider last 20%
% 
% if length(Amp)>100
% %     Amplitude = Amp(end);    
%     Amplitude = mean(Amp(AmpRange:end));    
% %     Period = t_peaks(end)-t_peaks(end-1);    
%     Period = mean(t_peaks(tpkRange+1:end)-t_peaks(tpkRange:end-1)); 
% else
%     Amplitude = 0;
%     Period = 0;
% end
% 
% Amplitude_store(k) = Amplitude;
% Period_store(k) = Period;
% 
% disp(k)

end


%% Percentage of synthetic protein per proteome

Hsum = HA+HB+HC;

total_aa = (R+TL_all)*n_R + C*n_C + Q*n_Q + Hsum*300;

R_frac = (R+TL_all)*n_R ./ total_aa;
C_frac = C*n_C ./ total_aa;
Q_frac = Q*n_Q ./ total_aa;
Hsum_frac = Hsum*300 ./ total_aa;

R_aa = (R+TL_all)*n_R;
C_aa = C*n_C;
Q_aa = Q*n_Q;
Hsum_aa = Hsum*300;


%% Plot individual population

% close
% xx=8;
% set(gcf, 'Position',  [900, 200, 550, 500])
% 
% plot(t, subpop(:,xx), 'LineWidth', 2)
% hold on
% 
% xlim([0 80])
% xlabel('Time / h')
% ylabel('Subpop')
% % ylim([0 1e9])
% 
% axis square


%% Plot single panel

% close
% set(0,'DefaultLineLineWidth',2.5);
% set(gcf, 'Position',  [800, 50, 600, 400])
% 
% plot(t, HA(:,1))
% hold on
% plot(t, HB(:,1))
% plot(t, HC(:,1))
% xlim([0 22])
% xticks([0, 10, 20])
% ylim([0 2500])
% yticks([0,1250,2500])
% xlabel('Time')
% ylabel({'Protein quantity'})
% grid on


%% State-specific proteins + GR/totalprotein

% xx=5;
% GRgreen = [26, 173, 27]/255;
% 
% close
% set(0,'DefaultLineLineWidth',2.5);
% set(gcf, 'Position',  [800, 50, 600, 700])
% 
% subplot(2,1,1)
% plot(t, HA(:,xx))
% hold on
% plot(t, HB(:,xx))
% plot(t, HC(:,xx))
% xlim([0 40])
% xticks([0, 10, 20])
% % ylim([0 9000])
% % yticks([0,4000,8000])
% ylabel({'Protein per cell'})
% grid on
% ax1=gca;
% 
% subplot(2,1,2)
% plot(t, HA(:,xx)+HB(:,xx)+HC(:,xx), 'k')
% hold on
% yyaxis right
% plot(t, GR(:,xx)./GR(:,end), 'color', GRgreen)
% ylim([0.2, 1])
% yticks([0.2,1])
% yticklabels({'0.2','1.0'})
% ylabel({'Relative','growth rate'})
% ax = gca;
% ax.YAxis(2).Color = GRgreen;
% yyaxis left
% xlim([0 25])
% xticks([0, 10, 20])
% ylim([0 9000])
% yticks([0,4000,8000])
% xlabel('Time/h')
% ylabel({'Protein per cell'})
% grid on
% ax2=gca;

% exportgraphics(ax1,'_Results/Repressilator_p4e3_3state_subpop14_UP.pdf','BackgroundColor','none','ContentType','vector')
% exportgraphics(ax2,'_Results/Repressilator_p4e3_3state_subpop14_DN.pdf','BackgroundColor','none','ContentType','vector')



%% Plot different proteins + their average for whole population - INTERMEDIATES

% close
% set(0,'DefaultLineLineWidth',1.5);
% set(0,'defaultAxesFontSize',24);
% 
% set(gcf, 'Position',  [500, 50, 1000, 700])
% 
% subplot(2,1,1)
% plot(t, HA_pop_percell)
% hold on
% plot(t, HB_pop_percell)
% plot(t, HC_pop_percell)
% xlim([0 120])
% ylim([0 9000])
% xticks([])
% yticks([0,4000,8000])
% ylabel('Protein per cell')
% grid on
% 
% yyaxis right
% a1 = area(t, subpop(:,1));
% a1.FaceColor = [1 1 1];
% a1.EdgeColor = [0 0 0];
% a1.LineWidth = 1;
% hold on
% 
% a8 = area(t,subpop(:,8));
% a8.FaceColor = E_yellow;
% a8.EdgeColor = [0 0 0];
% a8.LineWidth = 1;
% 
% a14 = area(t,subpop(:,14));
% a14.FaceColor = E_yellow;
% a14.EdgeColor = [0 0 0];
% a14.LineWidth = 1;
% 
% a27 = area(t,subpop(:,end));
% a27.FaceColor = M_red;
% a27.EdgeColor = [0 0 0];
% a27.LineWidth = 1;
% 
% alpha(0.1)
% 
% ylabel('Subpop. size')
% ax = gca;
% ax.YAxis(2).Color = 'k';
% ylim([0,1.1e9])
% yticks([0,1e9])
% yticklabels({'0','1e9'})
% yyaxis left
% ax1=gca;
% 
% subplot(2,1,2)
% plot(t, H_pop_percell, 'k')
% % hold on
% % yyaxis right
% % plot(t, rel_GR_pop_percell, 'color', GRgreen)
% % ylim([0.5, 1])
% % yticks([0.5,1])
% % yticklabels({'0.5','1.0'})
% % ylabel({'Relative','growth rate'})
% % ax = gca;
% % ax.YAxis(2).Color = GRgreen;
% % yyaxis left
% xlim([0 120])
% ylim([0 9000])
% xticks([0,20,40,60,80,100,120])
% yticks([0,4000,8000])
% xlabel('Time/h')
% ylabel('Protein per cell')
% grid on
% ax2=gca;

% exportgraphics(ax1,'_Results/Repressilator_p4e3_3state_pop_UP.png','Resolution',600) % SAVE
% exportgraphics(ax2,'_Results/Repressilator_p4e3_3state_pop_DN.png','Resolution',600) % SAVE



%% Paper plot - Population Protein

close
set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',24);

set(gcf, 'Position',  [500, 200, 507, 327])

plot(t, HA_pop_percell, 'Color', ABC_A)
hold on
plot(t, HB_pop_percell, 'Color', ABC_B)
plot(t, HC_pop_percell, 'Color', ABC_C)
% xlabel('Time/h')
xlim([0 25])
ylim([0 3e4])
xticks([0 25])

yyaxis right
a1 = area(t, subpop(:,1));
a1.FaceColor = [1 1 1];
a1.EdgeColor = [0 0 0];
a1.LineWidth = 1;
hold on

a27 = area(t,subpop(:,end));
a27.FaceColor = M_red;
a27.EdgeColor = [0 0 0];
a27.LineWidth = 1;

alpha(0.1)

ax = gca;
ax.YAxis(2).Color = 'k';
ylim([0,1.15e9])
yticks([0,1e9])
yyaxis left

% exportgraphics(gca,'_Results/PaperResults/rep_s3_pop_p1e5_1e4_z1e-3_5e-1.pdf','BackgroundColor','none','ContentType','vector')



%% Paper plot - GR and Energy per state

% close
% xx=2;
% GRgreen = [26, 173, 27]/255;
% set(0,'DefaultLineLineWidth',3);
% set(0,'defaultAxesFontSize',24);
% set(gcf, 'Position',  [500, 200, 442, 272])
% 
% % plot(t, Hsum_frac(:,xx), 'k')
% plot(t, e(:,xx), 'k')
% hold on
% 
% xlim([0 22])
% xticks([0 22])
% ylim([9 10.1])
% yticks([9 10])
% yline(mean(e(:,1)), 'k--', 'LineWidth', 2)
% 
% yyaxis right
% plot(t, GR(:,xx)./GR(:,end), 'color', GRgreen)
% ylim([0 0.4])
% yticks([0 0.4])
% yticklabels({'0','0.4'})
% ax = gca;
% ax.YAxis(2).Color = GRgreen;
% yyaxis left

% exportgraphics(gca,'_Results/PaperResults/rep_s3_subpop2_p1e5_1e4_z1e-3_5e-1.pdf','BackgroundColor','none','ContentType','vector')



%% SI plot - cells per state

% close
% f = figure;
% set(gcf, 'Position',  [700, 150, 800, 700])
% set(0,'DefaultLineLineWidth',7);
% set(0,'defaultAxesFontSize',16);
% 
% fig=1;
% 
% for i = 19:27
%     
% subplot(3, 3, fig)
% plot(t,subpop(:,i))
% 
% xlim([0 22])
% ylim([0 1e9])
% xlabel('')
% ylabel('')
% xticks([])
% yticks([])
% hold on
% axis square
% 
% fig = fig+1;
% 
% end



%% SI plot - oscillations

% close
% f = figure;
% set(gcf, 'Position',  [700, 150, 800, 700])
% set(0,'DefaultLineLineWidth',3);
% set(0,'defaultAxesFontSize',16);
% 
% fig=1;
% 
% for i = 19:27
%     
% subplot(3, 3, fig)
% plot(t,HA(:,i), 'Color', ABC_A)
% hold on
% plot(t,HB(:,i), 'Color', ABC_B)
% plot(t,HC(:,i), 'Color', ABC_C)
% 
% xlim([40 50])
% ylim([0 3e4])
% xlabel('')
% ylabel('')
% xticks([])
% yticks([])
% hold on
% axis square
% 
% fig = fig+1;
% 
% end


%% Plot certain subpop trajectories

% close
% set(gcf, 'Position',  [1000, 100, 900, 700])
% 
% plot(t,subpop(:,1))
% hold on
% plot(t,subpop(:,2))
% plot(t,subpop(:,5))
% plot(t,subpop(:,14))
% % plot(t,subpop(:,15))
% % plot(t,subpop(:,18))
% plot(t,subpop(:,27))
% 
% xlim([0 600])
% ylim([0 1e9])


%% Compare with isolated repressilator - simple protein plot

% close
% set(gcf, 'Position',  [600, 100, 800, 800])
% 
% plot(t, HA(:,1))
% hold on
% plot(t, HB(:,1))
% plot(t, HC(:,1))
% % xlim([0 30])
% % ylim([0 2e4])
% xlabel('Time / h')
% ylabel('Protein quantity')
% grid on
% axis square



%% Compare with isolated repressilator - Amplitude and Period

% set(0,'DefaultLineLineWidth',3);
% set(0,'defaultAxesFontSize',28);
% 
% Amplitude_store_EnzNut = [0,0,503.831800351622,724.419487834321,1017.55422043248,1410.02347160201,1934.94910262499,2632.52185980116,3549.31411765107,4734.30678901387,6228.93198081406,8048.63362866213,10155.5891502208,12431.4123732303,14677.7830928394,16672.4563199493,18255.3764821414,19362.8503703222,19990.6966910506,20150.7846689866,19857.5235265781,19132.2903449713,18011.2115955174,16551.4765287485,14834.5004568337,12964.5107393401,11057.2834652785,9216.33440610838,7502.45116934032,5926.62692369878];
% Period_store_EnzNut = [0,0,168.738582701330,174.218260122161,180.282180823340,186.679183252942,193.118978033170,199.272799005019,204.771299051400,209.196413836509,212.080195442779,212.948901664243,211.410809767742,207.350111941775,201.140503636636,193.702932146057,186.210568993812,179.505022900658,173.875829771312,169.231506927065,165.352120778216,162.027847017086,159.087662277625,156.376392234633,153.751300063505,151.074219955758,148.210761357275,145.034976703467,141.369389594562,136.920261965263];
% 
% Amplitude_store_MutBase = [0,0,503.851375500730,724.444182785651,1017.59974201327,1410.10576248965,1935.02524959548,2632.43018366788,3549.29511313085,4734.49580640116,6229.15401413620,8048.57789132695,10155.5432131031,12431.4058745770,14677.8285381190,16672.3906890165,18255.2124895421,19362.6003727297,19990.3013698522,20150.4062775510,19857.2734445137,19131.9571592144,18010.9460141893,16551.1793169173,14834.2127820103,12964.2386527431,11057.0769937876,9216.15243224185,7502.32512747111,5926.52858227003];
% Period_store_MutBase = [0,0,2.81232284357080,2.90364916902922,3.00471812331956,3.11134258661904,3.21866571762717,3.32123929916657,3.41288600141332,3.48662255038928,3.53467518953414,3.54914484045176,3.52352348177434,3.45587405049682,3.35232546725226,3.22842509539094,3.10350731460737,2.99173837081095,2.89793686727496,2.82052567471341,2.75586619898213,2.70047269007637,2.65144001446997,2.60626076021348,2.56252648171670,2.51790755447901,2.47018468202403,2.41722157686595,2.35615790178373,2.28202876389052];
% 
% close
% set(gcf, 'Position',  [800, 75, 750, 900])
% 
% subplot(2,1,1)
% 
% semilogx(induc_vec/60, Amplitude_store_EnzNut, '-', 'color', [222, 22, 22]/255)
% hold on
% semilogx(induc_vec/60, Amplitude_store_MutBase, 'kx')
% xlim([1e0, 1e4])
% xticks([1e0, 1e1, 1e2, 1e3, 1e4])
% ylim([0 2.4e4])
% yticks([0, 0.8e4 1.6e4, 2.4e4])
% ylabel('Amplitude')
% grid on
% 
% subplot(2,1,2)
% 
% semilogx(induc_vec/60, Period_store_EnzNut, '-', 'color', [222, 22, 22]/255)
% hold on
% semilogx(induc_vec/60, Period_store_MutBase*60, 'kx')
% xlim([1e0, 1e4])
% xticks([1e0, 1e1, 1e2, 1e3, 1e4])
% ylim([0 240])
% yticks([0, 80, 160, 240])
% xlabel('Induction Strength')
% ylabel('Period/h')
% grid on



%% Thesis extra plot(s)

% close
% set(gcf, 'Position',  [900, 200, 1100, 700])
% set(0,'DefaultLineLineWidth',3);
% set(0,'defaultAxesFontSize',22);
% 
% p_I_strength = [linspace(90,10,9),5];
% t_H1 = [292.9, 304, 327, 362, 413.4, 478.8, 605.4, 821.9, 1182, 1559.1];
% Harea = [5.85e5, 5.86e5, 5.90e5, 5.94e5, 6.00e5, 6.06e5, 6.12e5, 6.20e5, 6.24e5, 6.24e5];
% 
% subplot(1,2,1)
% 
% plot(p_I_strength, t_H1);
% yline(289.2)
% 
% xlim
% ylim([0 1500])
% yticks([0,500,1000,1500])
% xlabel({'Intermediate state','promoter strength'})
% ylabel('Time before crash/h')
% axis square
% grid on
% 
% subplot(1,2,2)
% 
% plot(p_I_strength, Harea);
% 
% ylim([5.7e5,6.3e5])
% yticks([5.7e5,6e5,6.3e5])
% % yticklabels({'5.5','6.5'})
% xlabel({'Intermediate state','promoter strength'})
% ylabel('Protein yield')
% axis square
% grid on












