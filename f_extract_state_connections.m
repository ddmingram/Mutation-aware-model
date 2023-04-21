%%% Additional script that converts information from 'state_connections'
%%% into more useable variables. The ODE for each state involves terms
%%% regarding (i) number of cells in connected upstream states, (ii)
%%% mutation probabilities from connected upstream states, (iii) mutation
%%% probabilities to downstream states. The cell array 'state_connections'
%%% has this info in an over-extensive and coded format. Hence,
%%% f_extract_state_connections transfers the required info into defined
%%% variables and converts 'coordinates' to numbered indices and
%%% 'z_matrix index pairs' into probabilities. For more contextual info,
%%% see "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

function [states_up,...
          z_values_up,...
          z_values_dn] = f_extract_state_connections(state_connections, active_z_matrix, s, d, n)
      
%%% Pre-allocate variables
[states_up, z_values_up, states_dn, z_values_dn] = deal(cell(n,1));

%%% Find upstream/downstream. For each state, store indicies of their up/dn
%%% states and z-values. Info for 'cells_dn' is not needed, but can be
%%% stored regardless.

for i=1:n

    num_up = size(state_connections{i,3}, 1);
    num_dn = size(state_connections{i,5}, 1);
    
    %%% Upstream
    for j = 1:num_up
        
        % States
        coord_up = state_connections{i,3}(j,:); % Get coord of first upstream state
        idx = f_coord_converter(coord_up, s, d); % What idx is that state?
        states_up{i} = [states_up{i}, idx];
        
        % z: values from state_connections act as indicies for z_mat
        state = state_connections{i,4}(j,1);
        dim = state_connections{i,4}(j,2);        
        z_values_up{i} = [z_values_up{i}, active_z_matrix(state,dim)];
    end    
    
    %%% Downstream
    for j = 1:num_dn
        
        % " "
        coord_dn = state_connections{i,5}(j,:);
        idx = f_coord_converter(coord_dn, s, d);
        states_dn{i} = [states_dn{i}, idx];
        
        % " "
        state = state_connections{i,6}(j,1);
        dim = state_connections{i,6}(j,2);        
        z_values_dn{i} = [z_values_dn{i}, active_z_matrix(state,dim)];
    end    
end