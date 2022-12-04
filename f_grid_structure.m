%%% Take in values of s and d...
%%% Return a full compilation of subpopulations with their coords,
%%% upstream/downstream subpops, and z_values.

function MyPop = f_grid_structure(s,d)

n = s^d;

%%
%%%%%%%%%%%%%%%%%%%%%
% Index cells for each dimension
%%%%%%%%%%%%%%%%%%%%%

coord_system = zeros(n,d);

for i = 1:d
    % 1. We want to cycle coords every s^(i-1)
    % ...eg when s=3,d=3, cycle every 3^2=9
    % 2. We then repeat each cycle s^(d-i) times
    % ...eg for s=3,d=4, when d=2 we repeat s^(4-2)=9 times 
    
    % Repeat elements that cycle every x --> repelem(0:s-1, x)
    % Repeat ABOVE y times --> repmat(repelem(0:s-1, x), [1,y])
    
    % Order coords 'backward's: x:-1:y
    
    coord = repmat( repelem(s-1:-1:0, s^(i-1)), [1,s^(d-i)] )';
    coord_system(:,i) = coord;
end



%%
%%%%%%%%%%%%%%%%%%%%%
% Create MyPop array
%%%%%%%%%%%%%%%%%%%%%

% Rows - 
% index | cell coord
% coords of uptream cells | z for upstream cells
% coords of downstream cells | z for downstream cells

MyPop = cell(n,6); % Pre-allocate
for i = 1:n
    MyPop(i,1) = {i};
    MyPop(i,2) = {coord_system(i,:)};
end



%%
%%%%%%%%%%%%%%%%%%%%%
% Create transitions FROM and TO
%%%%%%%%%%%%%%%%%%%%%

% 'partial coord' (pCoord) is one number within the full coordinate.
% For interest, CurrentState = s - pCoord

for i = 1:n         % For every subpop...
    for j = 1:d     % Check transitions for each dim   
        
        old_pCoord = MyPop{i,2}(j);             % For recording z_type
        
        % Find all UPSTREAM cells        
        new_pCoord = MyPop{i,2}(j);             % Take coordinate of current subpop + record pCoord of current dim
        
        while new_pCoord ~= s-1                 % Are backwards transitions possible?  
            % For next subpop, increment state
            new_pCoord = new_pCoord + 1;        % Increment
            pre_subpop = MyPop{i,2};            % Initialise coord for upstream subpop
            pre_subpop(j) = new_pCoord;         % Set coord
            MyPop{i,3} = [MyPop{i,3}; pre_subpop];  % Add to MyPop array   
            
            % For upstream, z_type relates to *current* subpop
            z_type = [old_pCoord+1, j]; % Rows = decreasing severity. Idx = pCoord+1
            MyPop{i,4} = [MyPop{i,4}; z_type];
        end
        
        % Find all DOWNSTREAM cells
        new_pCoord = MyPop{i,2}(j);
        
        while new_pCoord ~= 0
            % For next subpop, reduce state
            new_pCoord = new_pCoord - 1;        % Reduce state
            post_subpop = MyPop{i,2};
            post_subpop(j) = new_pCoord;
            MyPop{i,5} = [MyPop{i,5}; post_subpop];
            
            % For downstream, z_type relates to *new* subpop
            z_type = [new_pCoord+1, j];
            MyPop{i,6} = [MyPop{i,6}; z_type];
        end  
    end    
end

% disp(coord_system)