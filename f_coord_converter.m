%%% Script that takes a coordinate (reprenting a position in multidimensional
%%% space with a given d and s) and spits out the correct index number. Will
%%% likely work similar to a binary system.

%%%%%%%%%%%%%%%%%%%%%
% Counting where [0,0,...] = 1
%%%%%%%%%%%%%%%%%%%%%

% function index = a_coord_converter(coord,d,s)
% 
%     index = 1;
%     
%     for i = 1:d
%         unit_value = s^(i-1) * coord(i);
%         index = index + unit_value;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%
% Counting where [s-1,s-1,...] = 1
%%%%%%%%%%%%%%%%%%%%%

function index = f_coord_converter(coord,d,s)

    index = 1;
    
    for i = 1:d
        unit_value = s^(i-1) * (s-coord(i)-1);
        index = index + unit_value;
    end
end