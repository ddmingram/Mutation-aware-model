function y_act = convert_state_vec(y, bioreactor)

    cell_type_num = length(bioreactor.y_order);

    offset=0;
    
    for k=1:cell_type_num

        cell_type = bioreactor.y_order{k};

        if isfield(bioreactor,cell_type)
            state_num = size(fieldnames(bioreactor.(cell_type)),1);
            y_act.(cell_type) = y(1+offset:state_num+offset);
            offset = offset + state_num;

        else
            error('no such field: %s',cell_type)

        end
    end

end