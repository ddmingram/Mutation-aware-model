%%% Each subpop requires a unique set of part strengths, depending on the
%%% varied combination of dimensions and states we choose. This function
%%% takes the requied "d" and "s" values, and returns the part parameters
%%% with values for each subpop.

function [z_mat, parts_mat,...
            prom_plusA_vec,...
            prom_minusA_vec,...
            RBS_plusA_vec,...
            RBS_minusA_vec,...
            CDSA_vec,...
            prom_plusB_vec,...
            prom_minusB_vec,...
            RBS_plusB_vec,...
            RBS_minusB_vec,...
            CDSB_vec] = f_PartParameters_2(dim_vec, d, s, n,...
                                           parts_mat_all,...
                                           z_mat_all,...
                                           n_HA, n_HB)      
        
%% 1. Create vectors of active parts

%%% Part strengths
parts_mat = zeros(s,d);
count = 1;
for i=1:length(dim_vec)
    if dim_vec(i)==1
        if s==1 % Different indexing if all E cells
            parts_mat(1,count) = parts_mat_all(1,i);
        else
            parts_mat(1:s,count) = parts_mat_all([1:s-1,end], i);
            % E.g. if s=3, parts_mat_all indicies are 1,2,4
        end
        count=count+1;
    end
end

%%% Mutation probabilities
z_mat = zeros(s-1,d);
count = 1;
for i=1:length(dim_vec)
    if dim_vec(i)==1
        if s == 2
            z_mat(1,count) = z_mat_all(1,i);
            count=count+1;
        else
            z_mat(1:s-1,count) = z_mat_all([1:s-2,end], i);
            % E.g. if s=3, z_mat_all indicies are 1,3
            count=count+1;
        end
    end
end

% Adjust probabilities to fix total z across mutation events
% Important to compare effects of intermediates/part-specific!
z_mat = z_mat / sum(dim_vec);



%% 2. Allocate active parts (cols) to subpops in order (rows)

parts_combos = zeros(n,d);
for i = 1:d
    % For explanation, see grid_structure.m
    parts_combos(:,i) = repmat( repelem(parts_mat(1:end,i),...
                                   s^(i-1)), [s^(d-i),1] );
end



%% 3. Convert each column into correct part variable

ActiveDim=1; % For allocating active parts
x=1; % Part Type no.

%%%%%%%%%%%%%%% CONSTRUCT A %%%%%%%%%%%%%%%
if dim_vec(x) == 1 % Promoter
    prom_plusA_vec = parts_combos(:,ActiveDim)';
    prom_minusA_vec = prom_plusA_vec;
    ActiveDim=ActiveDim+1;
else
    prom_plusA_vec = [ones(n,1)*parts_mat_all(1,x)]';
    prom_minusA_vec = prom_plusA_vec;
end
x=x+1;

if dim_vec(x) == 1 % RBS
    RBS_plusA_vec = parts_combos(:,ActiveDim)';
    RBS_minusA_vec = RBS_plusA_vec;
    ActiveDim=ActiveDim+1;
else
    RBS_plusA_vec = [ones(n,1)*parts_mat_all(1,x)]';
    RBS_minusA_vec = RBS_plusA_vec;
end
x=x+1;

if dim_vec(x) == 1 % CDS
    CDSA_vec = parts_combos(:,ActiveDim)' .* n_HA;
    
    noHmutA_idx = []; % For H_avg, need subpops that mutate CDS
    for i=1:length(CDSA_vec)
        if CDSA_vec(i) == CDSA_vec(1) % If subpop doesn't contain mutated CDS
            noHmutA_idx = [noHmutA_idx,i]; % Index of subpop
        end
    end
else
    CDSA_vec = [ones(n,1) .* parts_mat_all(1,x) * n_HA]';
end
x=x+1;

%%%%%%%%%%%%%%% CONSTRUCT B %%%%%%%%%%%%%%%
if dim_vec(x) == 1 % Promoter
    prom_plusB_vec = parts_combos(:,ActiveDim)';
    prom_minusB_vec = prom_plusB_vec;
    ActiveDim=ActiveDim+1;
else
    prom_plusB_vec = [ones(n,1)*parts_mat_all(1,x)]';
    prom_minusB_vec = prom_plusB_vec;
end
x=x+1;

if dim_vec(x) == 1 % RBS
    RBS_plusB_vec = parts_combos(:,ActiveDim)';
    RBS_minusB_vec = RBS_plusB_vec;
    ActiveDim=ActiveDim+1;
else
    RBS_plusB_vec = [ones(n,1)*parts_mat_all(1,x)]';
    RBS_minusB_vec = RBS_plusB_vec;
end
x=x+1;

if dim_vec(x) == 1 % CDS
    CDSB_vec = parts_combos(:,ActiveDim)' .* n_HB;
    
    noHmutB_idx = []; % For H_avg, need subpops that mutate CDS
    for i=1:length(CDSB_vec)
        if CDSB_vec(i) == CDSB_vec(1)
            noHmutB_idx = [noHmutB_idx,i];
        end
    end
else
    CDSB_vec = [ones(n,1) .* parts_mat_all(1,x) * n_HB]';
end