%%% Additional script that labels each state systematically with a
%%% coordinate ('coord') of length 's', then creates the cell array
%%% state_connections that collects info on how each state is connected to
%%% one another. This is required for automatically assigning 'part' and
%%% 'z' values to each state in the correct order. This script works for
%%% any value of [s,d]. For more contextual info, see
%%% "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

function state_connections = f_allocate_coords(s,d,n)

%% Create systematic indices for each state

coord_system = zeros(n,d); % Initialise

for i = 1:d
    
    % Goal: create column of indices whose rows represent each state, one
    % dimension at a time. Order this top-down. E.g. for [s3,d2]: column
    % vectors would be:[2,1,0,2,1,0,2,1,0], [2,2,2,1,1,1,0,0,0].
    
    %%% 1. We want to cycle coords every s^(i-1)
    % -> x = repelem(s-1:-1:0, s^(i-1))
    % ...e.g. when s=3,d=3, cycle every 3^2=9
    
    %%% 2. We then repeat THAT cycle s^(d-i) times
    % repmat(x, [1,s^(d-i)])
    % ...e.g. for s=3,d=4: when d=2 we repeat s^(4-2)=9 times
    
    %%% 3. Order coords in reverse: a:-1:b
    
    coord = repmat( repelem(s-1:-1:0, s^(i-1)), [1,s^(d-i)] )';
    coord_system(:,i) = coord;
end



%% Create 'state_connections' cell array

%%% Rows: states in order. Columns: (1) index of state. (2) coord of state.
%%% (3) List of coords of upstream states. (4) List of index pairs that
%%% each specify a value from the z_matrix for the mutation probability of
%%% an upstream state. (5-6) Same as (3-4) but for downstream states.

state_connections = cell(n,6); % Initialise

%%% First two rows are simple
for i = 1:n
    state_connections(i,1) = {i}; % Row 1
    state_connections(i,2) = {coord_system(i,:)}; % Row 2
end

%%% 'From' (upstream) and 'To' (downstream coords are more involved.
%%% Here-in, 'partial coord' (pCoord) refers to a single digit within the
%%% full coord. CurrentState = s - pCoord

for i = 1:n % For every state
    
    for j = 1:d % Check transitions for each dimension
        
        old_pCoord = state_connections{i,2}(j); % For recording z_type
        
        %%% 1. Find all UPSTREAM states        
        new_pCoord = state_connections{i,2}(j); % Take coordinate of current subpop and record pCoord of current dim
        
        while new_pCoord ~= s-1 % Are backwards transitions possible?            
            new_pCoord = new_pCoord + 1; % For next subpop, increment state
            pre_subpop = state_connections{i,2}; % Initialise coord for upstream subpop
            pre_subpop(j) = new_pCoord; % Set coord
            state_connections{i,3} = [state_connections{i,3}; pre_subpop]; % Add to cell array   
            
            % For upstream, z_type relates to *current* subpop
            z_type = [old_pCoord+1, j]; % Rows = decreasing severity. Index = pCoord+1
            state_connections{i,4} = [state_connections{i,4}; z_type];
        end
        
        %%% 2. Find all DOWNSTREAM states
        new_pCoord = state_connections{i,2}(j);
        
        while new_pCoord ~= 0            
            new_pCoord = new_pCoord - 1; % For next subpop, reduce state
            post_subpop = state_connections{i,2};
            post_subpop(j) = new_pCoord;
            state_connections{i,5} = [state_connections{i,5}; post_subpop];
            
            % For downstream, z_type relates to *new* subpop
            z_type = [new_pCoord+1, j];
            state_connections{i,6} = [state_connections{i,6}; z_type];
        end
        
    end % Dimensions
    
end % States