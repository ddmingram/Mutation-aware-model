function dydt = RHS_U_CELL(t, y, p, bioreactor, buffer)
    
    s_conf_E = bioreactor.E_CELL;
    s_conf_U = bioreactor.U_CELL;        
    
    Ucell_idx = getStateIdx('Ucell', s_conf_U);     
    Ucell = y.U_CELL(Ucell_idx);    
    e_idx = getStateIdx('e', s_conf_U);
    e = y.U_CELL(e_idx);    
    TX_R_idx = getStateIdx('TX_R', s_conf_U);
    TX_R = y.U_CELL(TX_R_idx);    
    TX_C_idx = getStateIdx('TX_C', s_conf_U);
    TX_C = y.U_CELL(TX_C_idx);    
    TX_P_idx = getStateIdx('TX_P', s_conf_U);
    TX_P = y.U_CELL(TX_P_idx);    
    TX_Q_idx = getStateIdx('TX_Q', s_conf_U);
    TX_Q = y.U_CELL(TX_Q_idx);    
    m_R_idx = getStateIdx('m_R', s_conf_U);
    m_R = y.U_CELL(m_R_idx);    
    m_C_idx = getStateIdx('m_C', s_conf_U);
    m_C = y.U_CELL(m_C_idx);    
    m_P_idx = getStateIdx('m_P', s_conf_U);
    m_P = y.U_CELL(m_P_idx);    
    m_Q_idx = getStateIdx('m_Q', s_conf_U);
    m_Q = y.U_CELL(m_Q_idx);    
    TL_R_idx = getStateIdx('TL_R', s_conf_U);
    TL_R = y.U_CELL(TL_R_idx);    
    TL_C_idx = getStateIdx('TL_C', s_conf_U);
    TL_C = y.U_CELL(TL_C_idx);    
    TL_P_idx = getStateIdx('TL_P', s_conf_U);
    TL_P = y.U_CELL(TL_P_idx);    
    TL_Q_idx = getStateIdx('TL_Q', s_conf_U);
    TL_Q = y.U_CELL(TL_Q_idx);    
    R_idx = getStateIdx('R', s_conf_U);
    R = y.U_CELL(R_idx);    
    C_idx = getStateIdx('C', s_conf_U);
    C = y.U_CELL(C_idx);    
    P_idx = getStateIdx('P', s_conf_U);
    P = y.U_CELL(P_idx);    
    Q_idx = getStateIdx('Q', s_conf_U);
    Q = y.U_CELL(Q_idx);    
    TX_H_idx = getStateIdx('TX_H', s_conf_U);
    TX_H = y.U_CELL(TX_H_idx);    
    m_H_idx = getStateIdx('m_H', s_conf_U);
    m_H = y.U_CELL(m_H_idx);    
    TL_H_idx = getStateIdx('TL_H', s_conf_U);
    TL_H = y.U_CELL(TL_H_idx);    
    H_idx = getStateIdx('H', s_conf_U);
    H = y.U_CELL(H_idx);   
    
    %%%%%%%%%%%%%%%% CALCULATIONS
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
    GR_U = TL_rate * TL_all / p.mass;
    
    % Growth Rate E
    Ecell_idx = getStateIdx('Ecell', s_conf_E);
    Ecell = y.E_CELL(Ecell_idx);
    e_E_idx= getStateIdx('e', s_conf_E);
    e_E = y.E_CELL(e_E_idx);    
    TL_R_E_idx= getStateIdx('TL_R', s_conf_E);
    TL_R_E = y.E_CELL(TL_R_E_idx);
    TL_C_E_idx= getStateIdx('TL_C', s_conf_E);
    TL_C_E = y.E_CELL(TL_C_E_idx);
    TL_P_E_idx= getStateIdx('TL_P', s_conf_E);
    TL_P_E = y.E_CELL(TL_P_E_idx);
    TL_Q_E_idx= getStateIdx('TL_Q', s_conf_E);
    TL_Q_E = y.E_CELL(TL_Q_E_idx);
    TL_H_E_idx= getStateIdx('TL_H', s_conf_E);
    TL_H_E = y.E_CELL(TL_H_E_idx);
    
    TL_rate_E = (p.v_TL * e_E) / (p.K_TL + e_E);
    TL_all_E = TL_R_E + TL_C_E + TL_P_E + TL_Q_E + TL_H_E;    
    GR_E = TL_rate_E * TL_all_E / p.mass;

    %%%%%%%%%%%%%%%% ODEs
    
    dydt(length(y.U_CELL),1) = 0; % Setup ODE structure

    % Ucell    
    dydt(Ucell_idx) = Ecell*(p.zu)*GR_E + Ucell*(1-p.z)*GR_U - Ucell*(p.z)*GR_U + (Ucell * buffer);

    % e    
    dydt(e_idx) = epsilon - (gamma_all) - GR_U*e;

    % TX_R    
    dydt(TX_R_idx) = p.kb_TX*P - p.ku_TX*TX_R - w_R - GR_U*TX_R;
    % TX_C    
    dydt(TX_C_idx) = p.kb_TX*P - p.ku_TX*TX_C - w_C - GR_U*TX_C;
    % TX_P    
    dydt(TX_P_idx) = p.kb_TX*P - p.ku_TX*TX_P - w_P - GR_U*TX_P;
    % TX_Q    
    dydt(TX_Q_idx) = p.kb_TX*P - p.ku_TX*TX_Q - w_Q - GR_U*TX_Q;

    % m_R    
    dydt(m_R_idx) = w_R - p.kb_TL*m_R*R + p.ku_TL*TL_R + gamma_R - (GR_U+p.m_deg)*m_R;
    % m_C    
    dydt(m_C_idx) = w_C - p.kb_TL*m_C*R + p.ku_TL*TL_C + gamma_C - (GR_U+p.m_deg)*m_C;
    % m_P    
    dydt(m_P_idx) = w_P - p.kb_TL*m_P*R + p.ku_TL*TL_P + gamma_P - (GR_U+p.m_deg)*m_P;
    % m_Q    
    dydt(m_Q_idx) = w_Q - p.kb_TL*m_Q*R + p.ku_TL*TL_Q + gamma_Q - (GR_U+p.m_deg)*m_Q;

    % TL_R    
    dydt(TL_R_idx) = p.kb_TL*m_R*R - p.ku_TL*TL_R - gamma_R - GR_U*TL_R;
    % TL_C    
    dydt(TL_C_idx) = p.kb_TL*m_C*R - p.ku_TL*TL_C - gamma_C - GR_U*TL_C;
    % TL_P    
    dydt(TL_P_idx) = p.kb_TL*m_P*R - p.ku_TL*TL_P - gamma_P - GR_U*TL_P;
    % TL_Q    
    dydt(TL_Q_idx) = p.kb_TL*m_Q*R - p.ku_TL*TL_Q - gamma_Q - GR_U*TL_Q;

    % R    
    dydt(R_idx) = gamma_R + p.ku_TL*(TL_all-TL_H) - p.kb_TL*(m_all-m_H)*R + gamma_all + p.Ufactor*(p.RBS_minus*TL_H - p.RBS_plus*m_H*R) - GR_U*R;
    % C    
    dydt(C_idx) = gamma_C - GR_U*C;
    % P    
    dydt(P_idx) = gamma_P + p.ku_TX*(TX_all-TX_H) - 4*p.kb_TX*P + w_all + p.Ufactor*(p.prom_minus*TX_H - p.prom_plus*p.plasmid*P) - GR_U*P;
    % Q    
    dydt(Q_idx) = gamma_Q - GR_U*Q;

    % TX_H    
    dydt(TX_H_idx) = p.Ufactor*(p.prom_plus*p.plasmid*P - p.prom_minus*TX_H) - w_H - GR_U*TX_H;
    % m_H    
    dydt(m_H_idx) = w_H - p.Ufactor*(p.RBS_plus*m_H*R + p.RBS_minus*TL_H) + gamma_H - (GR_U+p.m_deg)*m_H;
    % TL_H    
    dydt(TL_H_idx) = p.Ufactor*(p.RBS_plus*m_H*R - p.RBS_minus*TL_H) - gamma_H - GR_U*TL_H;
    % H    
    dydt(H_idx) = gamma_H - GR_U*H;    
    
end