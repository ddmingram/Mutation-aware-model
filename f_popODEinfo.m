%%% Each subpop requires a different set of ODE variables, as it has a
%%% different set of upstream/downstream cells affecting it. This function
%%% identifies those cells and their linked probabilities

function [subpops_up,...
          z_values_up,...
          subpops_dn,...
          z_values_dn]...
          = f_popODEinfo(z_mat, MyPop, s, d, n) 
      
      
      
%% Pre-allocate

[subpops_up, z_values_up, subpops_dn, z_values_dn] = deal(cell(n,1));



%% Find upstream/downstream

% For each subpop, store indicies of their up/dn subpops + z_values
% I don't need cells_dn, but store anyway...

for i=1:n

    num_up = size(MyPop{i,3}, 1);
    num_dn = size(MyPop{i,5}, 1);
    
    % Upstream
    for j = 1:num_up
        %%% Subpops
        coord_up = MyPop{i,3}(j,:);             % Get coord of first upstream subpop
        idx = f_coord_converter(coord_up,d,s);    % What idx is that subpop?
        subpops_up{i} = [subpops_up{i}, idx];
        %%% z: values from MyPop act as indicies for z_mat
        state = MyPop{i,4}(j,1);
        dim = MyPop{i,4}(j,2);        
        z_values_up{i} = [z_values_up{i}, z_mat(state,dim)];
    end    
    
    % Downstream
    for j = 1:num_dn
        %%% Subpops
        coord_dn = MyPop{i,5}(j,:);
        idx = f_coord_converter(coord_dn,d,s);
        subpops_dn{i} = [subpops_dn{i}, idx];
        %%% z: values from MyPop act as indicies for z_mat
        state = MyPop{i,6}(j,1);
        dim = MyPop{i,6}(j,2);        
        z_values_dn{i} = [z_values_dn{i}, z_mat(state,dim)];
    end    
end