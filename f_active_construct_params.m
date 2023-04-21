%%% Additional script that obtains the 'active' parameters for part
%%% strengths and z_values, and distributes them across the states in the
%%% correct order. For more contextual info, see
%%% "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

function [active_z_matrix,...
          alpha_A_vec,...
          beta_A_vec,...
          alpha_B_vec,...
          beta_B_vec,...
          alpha_C_vec,...
          beta_C_vec] = f_active_construct_params(active_dims, s, d, n,...
                                                  part_matrix,...
                                                  z_matrix)     
                                       
%% 1. Create arrays of values associated with the active dimensions

%%% Part strengths
active_part_matrix = zeros(s,d);
count = 1;
for i = 1:length(active_dims)
    if active_dims(i) == 1
        if s==1 % Different indexing if all E cells
            active_part_matrix(1,count) = part_matrix(1,i);
        else
            active_part_matrix(1:s,count) = part_matrix([1:s-1,end], i);
            % E.g. if s=3, parts_mat_all indicies are 1,2,4
        end
        count=count+1;
    end
end

%%% Mutation probabilities
active_z_matrix = zeros(s-1,d);
count = 1;
for i = 1:length(active_dims)
    if active_dims(i) == 1
        if s == 2
            active_z_matrix(1,count) = z_matrix(1,i);
            count=count+1;
        else
            active_z_matrix(1:s-1,count) = z_matrix([1:s-2,end], i);
            % E.g. if s=3, z_mat_all indicies are 1,3
            count=count+1;
        end
    end
end


%% 2. Allocate active parts (cols) to subpops in order (rows)

parts_combos = zeros(n,d);

for i = 1:d
    % For explanation, see f_allocate_coords.m
    parts_combos(:,i) = repmat( repelem(active_part_matrix(1:end,i),...
                                   s^(i-1)), [s^(d-i),1] );
end


%% 3. Convert each column into correct part variable

%%% This section assigns the correct values to each active part manually
%%% and in order of the 'active_dims' vector [gene_A_prom, gene_A_RBS,
%%% gene_B_prom etc.]. It has code to account for up to three seperate
%%% genes. For >3 genes, more code is required in the same format as below.

ActiveDim = 1; % For allocating active parts
x = 1; % Part Type no.
num_genes = length(active_dims)/2; % '2': we consider prom/RBS for each gene

%%%%%%%%%%%%%%% CONSTRUCT A %%%%%%%%%%%%%%%
if active_dims(x) == 1 % Promoter
    alpha_A_vec = parts_combos(:,ActiveDim)';
    ActiveDim=ActiveDim+1;
else
    alpha_A_vec = [ones(n,1)*part_matrix(1,x)]';
end
x=x+1;

if active_dims(x) == 1 % RBS
    beta_A_vec = parts_combos(:,ActiveDim)';
    ActiveDim=ActiveDim+1;
else
    beta_A_vec = [ones(n,1)*part_matrix(1,x)]';
end

if num_genes > 1    
x=x+1;

%%%%%%%%%%%%%%% CONSTRUCT B %%%%%%%%%%%%%%%
if active_dims(x) == 1 % Promoter
    alpha_B_vec = parts_combos(:,ActiveDim)';
    ActiveDim=ActiveDim+1;
else
    alpha_B_vec = [ones(n,1)*part_matrix(1,x)]';
end
x=x+1;

if active_dims(x) == 1 % RBS
    beta_B_vec = parts_combos(:,ActiveDim)';
    ActiveDim=ActiveDim+1;
else
    beta_B_vec = [ones(n,1)*part_matrix(1,x)]';
end

if num_genes > 2    
x=x+1;

%%%%%%%%%%%%%%% CONSTRUCT C %%%%%%%%%%%%%%%
if active_dims(x) == 1 % Promoter
    alpha_C_vec = parts_combos(:,ActiveDim)';
    ActiveDim=ActiveDim+1;
else
    alpha_C_vec = [ones(n,1)*part_matrix(1,x)]';
end
x=x+1;

if active_dims(x) == 1 % RBS
    beta_C_vec = parts_combos(:,ActiveDim)';
else
    beta_C_vec = [ones(n,1)*part_matrix(1,x)]';
end

end % If >1 gene
end % If >2 genes