function idx = getStateIdx(name, state_config)

    if isfield(state_config, name)   
        f_names = fieldnames(state_config);  
        idx = find(strcmp(f_names,name));

    else
        error('unknown state name: %s', name)
    end
    
end
