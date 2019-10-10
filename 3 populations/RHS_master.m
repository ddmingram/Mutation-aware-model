% RHS_master

function dydt = RHS_master(t, y, p, bioreactor)

    y_act = convert_state_vec(y, bioreactor);
    
    Ecell_idx = getStateIdx('Ecell', bioreactor.E_CELL);
    Ecell = y_act.E_CELL(Ecell_idx);
    Ucell_idx = getStateIdx('Ucell', bioreactor.U_CELL);
    Ucell = y_act.U_CELL(Ucell_idx);
    Mcell_idx = getStateIdx('Mcell', bioreactor.M_CELL);
    Mcell = y_act.M_CELL(Mcell_idx);
    
    buffer = p.N - (Ecell + Ucell + Mcell);

    dydt_E_CELL = RHS_E_CELL(t, y_act, p, bioreactor, buffer);
    dydt_U_CELL = RHS_U_CELL(t, y_act, p, bioreactor, buffer); 
    dydt_M_CELL = RHS_M_CELL(t, y_act, p, bioreactor, buffer);    
        
    dydt = [dydt_E_CELL; dydt_U_CELL; dydt_M_CELL];
    
end






