% RHS_master

function dydt = RHS_master(t, y, p, bioreactor)

    y_act = convert_state_vec(y, bioreactor);
    
    Ecell_idx = getStateIdx('Ecell', bioreactor.E_CELL);
    Ecell = y_act.E_CELL(Ecell_idx);
    Mpcell_idx = getStateIdx('Mpcell', bioreactor.MP_CELL);
    Mpcell = y_act.MP_CELL(Mpcell_idx);
    Mrcell_idx = getStateIdx('Mrcell', bioreactor.MR_CELL);
    Mrcell = y_act.MR_CELL(Mrcell_idx);
    Mprcell_idx = getStateIdx('Mprcell', bioreactor.MPR_CELL);
    Mprcell = y_act.MPR_CELL(Mprcell_idx);
    
    buffer = p.N - (Ecell + Mpcell + Mrcell + Mprcell);

    dydt_E_CELL = RHS_E_CELL(t, y_act, p, bioreactor, buffer);
    dydt_MP_CELL = RHS_MP_CELL(t, y_act, p, bioreactor, buffer);
    dydt_MR_CELL = RHS_MR_CELL(t, y_act, p, bioreactor, buffer);    
    dydt_MPR_CELL = RHS_MPR_CELL(t, y_act, p, bioreactor, buffer);    
        
    dydt = [dydt_E_CELL; dydt_MP_CELL; dydt_MR_CELL; dydt_MPR_CELL];
    
end






