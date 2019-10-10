function dydt = RHS_MPR_CELL(t, y, p, bioreactor, buffer)    

    s_conf_Mp = bioreactor.MP_CELL;
    s_conf_Mr = bioreactor.MR_CELL;
    s_conf_Mpr = bioreactor.MPR_CELL;
    
    %%%%%%%%%%%%%%%% CALCULATIONS
    Mprcell_idx = getStateIdx('Mprcell', s_conf_Mpr);     
    Mprcell = y.MPR_CELL(Mprcell_idx);    
    e_idx = getStateIdx('e', s_conf_Mpr);
    e = y.MPR_CELL(e_idx);    
    TX_R_idx = getStateIdx('TX_R', s_conf_Mpr);
    TX_R = y.MPR_CELL(TX_R_idx);    
    TX_C_idx = getStateIdx('TX_C', s_conf_Mpr);
    TX_C = y.MPR_CELL(TX_C_idx);    
    TX_P_idx = getStateIdx('TX_P', s_conf_Mpr);
    TX_P = y.MPR_CELL(TX_P_idx);    
    TX_Q_idx = getStateIdx('TX_Q', s_conf_Mpr);
    TX_Q = y.MPR_CELL(TX_Q_idx);    
    m_R_idx = getStateIdx('m_R', s_conf_Mpr);
    m_R = y.MPR_CELL(m_R_idx);    
    m_C_idx = getStateIdx('m_C', s_conf_Mpr);
    m_C = y.MPR_CELL(m_C_idx);    
    m_P_idx = getStateIdx('m_P', s_conf_Mpr);
    m_P = y.MPR_CELL(m_P_idx);    
    m_Q_idx = getStateIdx('m_Q', s_conf_Mpr);
    m_Q = y.MPR_CELL(m_Q_idx);    
    TL_R_idx = getStateIdx('TL_R', s_conf_Mpr);
    TL_R = y.MPR_CELL(TL_R_idx);    
    TL_C_idx = getStateIdx('TL_C', s_conf_Mpr);
    TL_C = y.MPR_CELL(TL_C_idx);    
    TL_P_idx = getStateIdx('TL_P', s_conf_Mpr);
    TL_P = y.MPR_CELL(TL_P_idx);    
    TL_Q_idx = getStateIdx('TL_Q', s_conf_Mpr);
    TL_Q = y.MPR_CELL(TL_Q_idx);    
    R_idx = getStateIdx('R', s_conf_Mpr);
    R = y.MPR_CELL(R_idx);    
    C_idx = getStateIdx('C', s_conf_Mpr);
    C = y.MPR_CELL(C_idx);    
    P_idx = getStateIdx('P', s_conf_Mpr);
    P = y.MPR_CELL(P_idx);    
    Q_idx = getStateIdx('Q', s_conf_Mpr);
    Q = y.MPR_CELL(Q_idx);    
    TX_H_idx = getStateIdx('TX_H', s_conf_Mpr);
    TX_H = y.MPR_CELL(TX_H_idx);    
    m_H_idx = getStateIdx('m_H', s_conf_Mpr);
    m_H = y.MPR_CELL(m_H_idx);    
    TL_H_idx = getStateIdx('TL_H', s_conf_Mpr);
    TL_H = y.MPR_CELL(TL_H_idx);    
    H_idx = getStateIdx('H', s_conf_Mpr);
    H = y.MPR_CELL(H_idx);   

    % Building blocks
    epsilon = C * p.v_e * p.s / (p.K_e + p.s);
    % Transcription
    TX_rate = (p.v_TX * e) / (p.K_TX + e);
    IQ = 1 / (1 + (Q/p.K_Q)^p.hQ);
    w_R = TX_R * TX_rate / p.n_R;
    w_C = TX_C * TX_rate / p.n_C;
    w_P = TX_P * TX_rate / p.n_P;
    w_Q = TX_Q * TX_rate / p.n_Q * IQ;
    w_H = TX_H * TX_rate / p.n_H;
    % Translation
    TL_rate = (p.v_TL * e) / (p.K_TL + e);
    gamma_R = TL_R * TL_rate / (p.n_R/3);
    gamma_C = TL_C * TL_rate / (p.n_C/3);
    gamma_P = TL_P * TL_rate / (p.n_P/3);
    gamma_Q = TL_Q * TL_rate / (p.n_Q/3);
    gamma_H = TL_H * TL_rate / (p.n_H/3);
    % Summations
    TX_all = TX_R + TX_C + TX_P + TX_Q + TX_H;
    w_all = w_R + w_C + w_P + w_Q + w_H;
    m_all = m_R + m_C + m_P + m_Q + m_H;
    TL_all = TL_R + TL_C + TL_P + TL_Q + TL_H;
    gamma_all = gamma_R + gamma_C + gamma_P + gamma_Q + gamma_H;
    % Growth Rate
    GR_MP = TL_rate * TL_all / p.mass;
    
    % Growth Rate Mp
    Mpcell_idx = getStateIdx('Mpcell', s_conf_Mp);
    Mpcell = y.MP_CELL(Mpcell_idx);
    e_Mp_idx= getStateIdx('e', s_conf_Mp);
    e_Mp = y.MP_CELL(e_Mp_idx);    
    TL_R_Mp_idx= getStateIdx('TL_R', s_conf_Mp);
    TL_R_Mp = y.MP_CELL(TL_R_Mp_idx);
    TL_C_Mp_idx= getStateIdx('TL_C', s_conf_Mp);
    TL_C_Mp = y.MP_CELL(TL_C_Mp_idx);
    TL_P_Mp_idx= getStateIdx('TL_P', s_conf_Mp);
    TL_P_Mp = y.MP_CELL(TL_P_Mp_idx);
    TL_Q_Mp_idx= getStateIdx('TL_Q', s_conf_Mp);
    TL_Q_Mp = y.MP_CELL(TL_Q_Mp_idx);
    TL_H_Mp_idx= getStateIdx('TL_H', s_conf_Mp);
    TL_H_Mp = y.MP_CELL(TL_H_Mp_idx);    
    
    TL_rate_Mp = (p.v_TL * e_Mp) / (p.K_TL + e_Mp);    
    TL_all_Mp = TL_R_Mp + TL_C_Mp + TL_P_Mp + TL_Q_Mp + TL_H_Mp;
    GR_Mp = TL_rate_Mp * TL_all_Mp / p.mass;
    
    % Growth Rate Mr    
    Mrcell_idx = getStateIdx('Mrcell', s_conf_Mr);
    Mrcell = y.MR_CELL(Mrcell_idx);
    e_Mr_idx= getStateIdx('e', s_conf_Mr);
    e_Mr = y.MR_CELL(e_Mr_idx);    
    TL_R_Mr_idx= getStateIdx('TL_R', s_conf_Mr);
    TL_R_Mr = y.MR_CELL(TL_R_Mr_idx);
    TL_C_Mr_idx= getStateIdx('TL_C', s_conf_Mr);
    TL_C_Mr = y.MR_CELL(TL_C_Mr_idx);
    TL_P_Mr_idx= getStateIdx('TL_P', s_conf_Mr);
    TL_P_Mr = y.MR_CELL(TL_P_Mr_idx);
    TL_Q_Mr_idx= getStateIdx('TL_Q', s_conf_Mr);
    TL_Q_Mr = y.MR_CELL(TL_Q_Mr_idx);
    TL_H_Mr_idx= getStateIdx('TL_H', s_conf_Mr);
    TL_H_Mr = y.MR_CELL(TL_H_Mr_idx);
    
    TL_rate_Mr = (p.v_TL * e_Mr) / (p.K_TL + e_Mr);
    TL_all_Mr = TL_R_Mr + TL_C_Mr + TL_P_Mr + TL_Q_Mr + TL_H_Mr;
    GR_Mr = TL_rate_Mr * TL_all_Mr / p.mass;
    
    %%%%%%%%%%%%%%%% ODEs
    
    dydt(length(y.MPR_CELL),1) = 0; % Setup ODE structure

    % Mprcell    
    dydt(Mprcell_idx) = Mpcell*p.zr*GR_Mp + Mrcell*p.zp*GR_Mr + Mprcell*(1)*GR_MP + (Mprcell * buffer);

    % e    
    dydt(e_idx) = epsilon - (gamma_all) - GR_MP*e;

    % TX_R    
    dydt(TX_R_idx) = p.kb_TX*P - p.ku_TX*TX_R - w_R - GR_MP*TX_R;
    % TX_C    
    dydt(TX_C_idx) = p.kb_TX*P - p.ku_TX*TX_C - w_C - GR_MP*TX_C;
    % TX_P    
    dydt(TX_P_idx) = p.kb_TX*P - p.ku_TX*TX_P - w_P - GR_MP*TX_P;
    % TX_Q    
    dydt(TX_Q_idx) = p.kb_TX*P - p.ku_TX*TX_Q - w_Q - GR_MP*TX_Q;

    % m_R    
    dydt(m_R_idx) = w_R - p.kb_TL*m_R*R + p.ku_TL*TL_R + gamma_R - (GR_MP+p.m_deg)*m_R;
    % m_C    
    dydt(m_C_idx) = w_C - p.kb_TL*m_C*R + p.ku_TL*TL_C + gamma_C - (GR_MP+p.m_deg)*m_C;
    % m_P    
    dydt(m_P_idx) = w_P - p.kb_TL*m_P*R + p.ku_TL*TL_P + gamma_P - (GR_MP+p.m_deg)*m_P;
    % m_Q    
    dydt(m_Q_idx) = w_Q - p.kb_TL*m_Q*R + p.ku_TL*TL_Q + gamma_Q - (GR_MP+p.m_deg)*m_Q;

    % TL_R    
    dydt(TL_R_idx) = p.kb_TL*m_R*R - p.ku_TL*TL_R - gamma_R - GR_MP*TL_R;
    % TL_C    
    dydt(TL_C_idx) = p.kb_TL*m_C*R - p.ku_TL*TL_C - gamma_C - GR_MP*TL_C;
    % TL_P    
    dydt(TL_P_idx) = p.kb_TL*m_P*R - p.ku_TL*TL_P - gamma_P - GR_MP*TL_P;
    % TL_Q    
    dydt(TL_Q_idx) = p.kb_TL*m_Q*R - p.ku_TL*TL_Q - gamma_Q - GR_MP*TL_Q;

    % R    
    dydt(R_idx) = gamma_R + p.ku_TL*(TL_all-TL_H) - p.kb_TL*(m_all-m_H)*R + gamma_all + p.Mfactor*(p.RBS_minus*TL_H - p.RBS_plus*m_H*R) - GR_MP*R;
    % C    
    dydt(C_idx) = gamma_C - GR_MP*C;
    % P    
    dydt(P_idx) = gamma_P + p.ku_TX*(TX_all-TX_H) - 4*p.kb_TX*P + w_all + p.Mfactor*(p.prom_minus*TX_H - p.prom_plus*p.plasmid*P) - GR_MP*P;
    % Q    
    dydt(Q_idx) = gamma_Q - GR_MP*Q;

    % TX_H    
    dydt(TX_H_idx) = p.Mfactor*(p.prom_plus*p.plasmid*P - p.prom_minus*TX_H) - w_H - GR_MP*TX_H;
    % m_H    
    dydt(m_H_idx) = w_H - p.Mfactor*(p.RBS_plus*m_H*R + p.RBS_minus*TL_H )+ gamma_H - (GR_MP+p.m_deg)*m_H;
    % TL_H    
    dydt(TL_H_idx) = p.Mfactor*(p.RBS_plus*m_H*R - p.RBS_minus*TL_H) - gamma_H - GR_MP*TL_H;
    % H    
    dydt(H_idx) = gamma_H - GR_MP*H;    
    
end