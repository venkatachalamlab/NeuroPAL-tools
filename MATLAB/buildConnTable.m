function [conn_table, class_table] = buildConnTable(conn_edges, conn_i, neurons, varargin)
%BUILDCONNTABLE Build a connectome table for Cytoscape.
%
%   Input:
%   conn_edges   - the complete connectome
%   conn_i       - the list of edges for the new connectome
%   neurons      - the list of neurons for the new connectome
%   ci_signs     - the confidence interval sign for the list of neurons
%                  + = excitation, - = inhibition
%   is_add_motor - should we add the motor neurons?

% Initialize the connectome info.
sources = unique(conn_edges.Source(conn_i));
targets = unique(conn_edges.Target(conn_i));

% Initialize the connectome subset.
% Keep: Source, Target, Weight, source_class, target_class, Type.
conn_table = conn_edges(conn_i,[1:5,end]);

%% Construct the connectome subset.
% Tranlsate the connections to the new table format.
chemical = categorical("chemical");
electrical = categorical("electrical");
premotor_inter = categorical("pre-motor inter");
premotor_classes = categorical({'AVA', 'AVB', 'AVD', 'AVE', 'PVC'});
conn_table.source_type(:) = categorical(nan); % sensory/inter/motor/poly
conn_table.target_type(:) = categorical(nan); % sensory/inter/motor/poly
conn_table.source_NT(:) = categorical(nan); % Ach/Glut/GABA or mix.
conn_table.edge_mNT(:) = categorical(nan); % Ach/Glut/GABA or mix.
for i = 1:length(conn_i)
    j = conn_i(i);
    
    % Determine the source type.
    if any(conn_edges.source_class(j) == premotor_classes)
        conn_table.source_type(i) = premotor_inter;
    elseif conn_edges.source_poly(j)
        conn_table.source_type(i) = categorical("polymodal");
    elseif conn_edges.source_sensory(j)
        conn_table.source_type(i) = categorical("sensory");
    elseif conn_edges.source_inter(j)
        conn_table.source_type(i) = categorical("inter");
    elseif conn_edges.source_motor(j)
        conn_table.source_type(i) = categorical("motor");
    end
    
    % Determine the target type.
    if any(conn_edges.target_class(j) == premotor_classes)
        conn_table.target_type(i) = premotor_inter;
    elseif conn_edges.target_poly(j)
        conn_table.target_type(i) = categorical("polymodal");
    elseif conn_edges.target_sensory(j)
        conn_table.target_type(i) = categorical("sensory");
    elseif conn_edges.target_inter(j)
        conn_table.target_type(i) = categorical("inter");
    elseif conn_edges.target_motor(j)
        conn_table.target_type(i) = categorical("motor");
    end
    
    % Determine the source classical neurotransmitter.
    NT_name = [];
    if  conn_edges.acetylcholine(j)
        NT_name = [NT_name 'Ach'];
    end
    if  conn_edges.glutamate(j)
        NT_name = [NT_name 'Glut'];
    end
    if  conn_edges.GABA(j)
        NT_name = [NT_name 'GABA'];
    end
    if ~isempty(NT_name)
        conn_table.source_NT(i) = categorical(cellstr(NT_name));
    end
    
    % Determine the classical metabotropic communication.
    if conn_edges.Type(j) == chemical
        mNT_name = [];
        if conn_edges.acetylcholine(j) && conn_edges.target_GAR(j)
            mNT_name = [mNT_name 'Gar'];
        end
        if conn_edges.glutamate(j) && conn_edges.target_MGL(j)
            mNT_name = [mNT_name 'Mgl'];
        end
        if conn_edges.GABA(j) && conn_edges.target_GBB(j)
            mNT_name = [mNT_name 'Gbb'];
        end
        if ~isempty(mNT_name)
            conn_table.edge_mNT(i) = categorical(cellstr(mNT_name));
        end
    end
end

% Add in any missing neurons.
% Note: some neurons may be missing because they're not connected to the
% primary network (of included neurons).
conn_neurons = union(sources, targets);
missing_sources = setdiff(conn_neurons, sources);
for i = 1:length(missing_sources)

    % Initialize the source row.
    source_i = find(conn_edges.Source == missing_sources(i), 1);
    source_row = conn_edges(source_i,:);
    
    % Initialize the source info.
    % Note: represent the missing neuron.
    j = size(conn_table,1) + 1;
    conn_table.Source(j) = source_row.Source;
    conn_table.source_class(j) = source_row.source_class;
    conn_table.edge_mNT(j) = categorical(nan);
    
    % Compute the chemical weights
    conn_table.Type(j) = categorical(nan);
    conn_table.Weight(j) = 0;
    
    % Initialize the target info.
    conn_table.Target(j) = categorical(nan);
    conn_table.target_class(j) = categorical(nan);
    conn_table.target_type(j) = categorical(nan);
    
    % Determine the source type.
    if any(source_row.source_class == premotor_classes)
        conn_table.source_type(j) = premotor_inter;
    elseif source_row.source_poly
        conn_table.source_type(j) = categorical("polymodal");
    elseif source_row.source_sensory
        conn_table.source_type(j) = categorical("sensory");
    elseif source_row.source_inter
        conn_table.source_type(j) = categorical("inter");
    elseif source_row.source_motor
        conn_table.source_type(j) = categorical("motor");
    end
    
    % Determine the classical neurotransmitter.
    NT_name = [];
    if source_row.acetylcholine
        NT_name = [NT_name 'Ach'];
    end
    if source_row.glutamate
        NT_name = [NT_name 'Glut'];
    end
    if source_row.GABA
        NT_name = [NT_name 'GABA'];
    end
    if ~isempty(NT_name)
        conn_table.source_NT(j) = categorical(cellstr(NT_name));
    end
end


%% Add the connectome's motor neuron targets.
if length(varargin) > 1 && varargin{2}
    
    % Initialize the motor neuron target group.
    motor = categorical("motor");
    motor_neurons = categorical("Motor Neurons");
        
    % Add the connectome's motor neurons.
    conn_sources = unique(conn_table.Source);
    for i = 1:length(conn_sources)
        
        %% Add the chemical motor synapses.
        % Add the chemical synapses TO motor neurons.
        chemical_i = find(conn_edges.Source == conn_sources(i) & ...
            conn_edges.target_motor & ...
            conn_edges.Type == chemical);
        if ~isempty(chemical_i)
            j = size(conn_table,1) + 1;
            
            % Initialize the source row.
            source_i = find(conn_table.Source == conn_sources(i), 1);
            source_row = conn_table(source_i,:);
            
            % Initialize the source info.
            conn_table.Source(j) = conn_sources(i);
            conn_table.source_class(j) = source_row.source_class;
            conn_table.source_type(j) = source_row.source_type;
            conn_table.source_NT(j) = source_row.source_NT;
            conn_table.edge_mNT(j) = categorical(nan);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_edges.Weight(chemical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = motor_neurons;
            conn_table.target_class(j) = motor_neurons;
            conn_table.target_type(j) = motor;
        end

        % Add the chemical synapses FROM motor neurons.
        chemical_i = find(conn_edges.Target == conn_sources(i) & ...
            conn_edges.source_motor & ...
            conn_edges.Type == chemical);
        if ~isempty(chemical_i)
            j = size(conn_table,1) + 1;
            
            % Initialize the source info.
            source_i = find(conn_table.Source == conn_sources(i), 1);
            source_row = conn_table(source_i,:);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_edges.Weight(chemical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = motor_neurons;
            conn_table.source_class(j) = motor_neurons;
            conn_table.source_type(j) = motor;
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
            
            % Initialize the target info.
            conn_table.Target(j) = conn_sources(i);
            conn_table.target_class(j) = source_row.source_class;
            conn_table.target_type(j) = source_row.source_type;
        end
        
        %% Add the electrical motor synapses.
        % Add the electrical synapses TO motor neurons.
        electrical_i = find(conn_edges.Source == conn_sources(i) & ...
            conn_edges.target_motor & ...
            conn_edges.Type == electrical);
        if ~isempty(electrical_i)
            j = size(conn_table,1) + 1;
            
            % Initialize the source info.
            source_i = find(conn_table.Source == conn_sources(i), 1);
            source_row = conn_table(source_i,:);
            
            % Initialize the source info.
            conn_table.Source(j) = conn_sources(i);
            conn_table.source_class(j) = source_row.source_class;
            conn_table.source_type(j) = source_row.source_type;
            conn_table.source_NT(j) = source_row.source_NT;
            conn_table.edge_mNT(j) = categorical(nan);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_edges.Weight(electrical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = motor_neurons;
            conn_table.target_class(j) = motor_neurons;
            conn_table.target_type(j) = motor;
        end
        
        % Add the electrical synapses FROM motor neurons.
        electrical_i = find(conn_edges.Target == conn_sources(i) & ...
            conn_edges.source_motor & ...
            conn_edges.Type == electrical);
        if ~isempty(electrical_i)
            j = size(conn_table,1) + 1;
            
            % Initialize the source info.
            source_i = find(conn_table.Source == conn_sources(i), 1);
            source_row = conn_table(source_i,:);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_edges.Weight(electrical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = motor_neurons;
            conn_table.source_class(j) = motor_neurons;
            conn_table.source_type(j) = motor;
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
            
            % Initialize the target info.
            conn_table.Target(j) = conn_sources(i);
            conn_table.target_class(j) = source_row.source_class;
            conn_table.target_type(j) = source_row.source_type;
        end
    end

    %% Combine the pre-motor neurons into a group.
    % Initialize the pre-motor neuron group.
    premotor_neurons = categorical("Pre-Motor Neurons");
    conn_neurons = union( ...
        conn_table.Source(conn_table.source_type ~= premotor_inter), ...
        conn_table.Target(conn_table.target_type ~= premotor_inter));
    for i = 1:length(conn_neurons)
        
        %% Add the chemical premotor synapses.
        % Add the chemical synapses TO premotor neurons.
        chemical_i = find(conn_table.Source == conn_neurons(i) & ...
            conn_table.target_type == premotor_inter & ...
            conn_table.Type == chemical);
        if ~isempty(chemical_i)
            
            % Initialize the source info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(chemical_i(1),:);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_table.Weight(chemical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = premotor_neurons;
            conn_table.target_class(j) = premotor_neurons;
            conn_table.target_type(j) = premotor_inter;
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        % Add the chemical synapses FROM premotor neurons.
        chemical_i = find(conn_table.Target == conn_neurons(i) & ...
            conn_table.source_type == premotor_inter & ...
            conn_table.Type == chemical);
        if ~isempty(chemical_i)
            
            % Initialize the target info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(chemical_i(1),:);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_table.Weight(chemical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = premotor_neurons;
            conn_table.source_class(j) = premotor_neurons;
            conn_table.source_type(j) = premotor_inter;
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        %% Add the electrical premotor synapses.
        % Add the electrical synapses TO premotor neurons.
        electrical_i = find(conn_table.Source == conn_neurons(i) & ...
            conn_table.target_type == premotor_inter & ...
            conn_table.Type == electrical);
        if ~isempty(electrical_i)
            
            % Initialize the source info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(electrical_i(1),:);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_table.Weight(electrical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = premotor_neurons;
            conn_table.target_class(j) = premotor_neurons;
            conn_table.target_type(j) = premotor_inter;
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        % Add the electrical synapses FROM premotor neurons.
        electrical_i = find(conn_table.Target == conn_neurons(i) & ...
            conn_table.source_type == premotor_inter & ...
            conn_table.Type == electrical);
        if ~isempty(electrical_i)
            
            % Initialize the target info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(electrical_i(1),:);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_table.Weight(electrical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = premotor_neurons;
            conn_table.source_class(j) = premotor_neurons;
            conn_table.source_type(j) = premotor_inter;
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
    end
    
    %% Combine the locomotion neurons into a group.
    % Initialize the locomotion neuron group.
    locomotion_neurons = categorical("Locomotion Neurons");
    conn_neurons = union( ...
        conn_table.Source(conn_table.source_type ~= premotor_inter & ...
        conn_table.source_type ~= motor), ...
        conn_table.Target(conn_table.target_type ~= premotor_inter & ...
        conn_table.target_type ~= motor));
    for i = 1:length(conn_neurons)
        
        %% Add the chemical locomotion synapses.
        % Add the chemical synapses TO locomotion neurons.
        chemical_i = find(conn_table.Source == conn_neurons(i) & ...
            (conn_table.Target == premotor_neurons | ...
            conn_table.Target == motor_neurons) & ...
            conn_table.Type == chemical);
        if ~isempty(chemical_i)
            
            % Initialize the source info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(chemical_i(1),:);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_table.Weight(chemical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = locomotion_neurons;
            conn_table.target_class(j) = locomotion_neurons;
            conn_table.target_type(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        % Add the chemical synapses FROM locomotion neurons.
        chemical_i = find(conn_table.Target == conn_neurons(i) & ...
            (conn_table.Source == premotor_neurons | ...
            conn_table.Source == motor_neurons) & ...
            conn_table.Type == chemical);
        if ~isempty(chemical_i)
            
            % Initialize the target info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(chemical_i(1),:);
            
            % Compute the chemical weights
            conn_table.Type(j) = chemical;
            conn_table.Weight(j) = sum(conn_table.Weight(chemical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = locomotion_neurons;
            conn_table.source_class(j) = locomotion_neurons;
            conn_table.source_type(j) = categorical(nan);
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        %% Add the electrical locomotion synapses.
        % Add the electrical synapses TO locomotion neurons.
        electrical_i = find(conn_table.Source == conn_neurons(i) & ...
            (conn_table.Target == premotor_neurons | ...
            conn_table.Target == motor_neurons) & ...
            conn_table.Type == electrical);
        if ~isempty(electrical_i)
            
            % Initialize the source info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(electrical_i(1),:);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_table.Weight(electrical_i));
            
            % Initialize the target info.
            conn_table.Target(j) = locomotion_neurons;
            conn_table.target_class(j) = locomotion_neurons;
            conn_table.target_type(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
        
        % Add the electrical synapses FROM locomotion neurons.
        electrical_i = find(conn_table.Target == conn_neurons(i) & ...
            (conn_table.Source == premotor_neurons | ...
            conn_table.Source == motor_neurons) & ...
            conn_table.Type == electrical);
        if ~isempty(electrical_i)
            
            % Initialize the target info.
            j = size(conn_table,1) + 1;
            conn_table(j,:) = conn_table(electrical_i(1),:);
            
            % Compute the electrical weights
            conn_table.Type(j) = electrical;
            conn_table.Weight(j) = sum(conn_table.Weight(electrical_i));
            
            % Initialize the source info.
            conn_table.Source(j) = locomotion_neurons;
            conn_table.source_class(j) = locomotion_neurons;
            conn_table.source_type(j) = categorical(nan);
            conn_table.source_NT(j) = categorical(nan);
            conn_table.edge_mNT(j) = categorical(nan);
        end
    end
end

% Add the confidence interval signs to the source neurons.
if ~isempty(varargin)
    
    % Initilialize the confidence interval neurons & signs.
    ci_signs = varargin{1};
    conn_table.source_sign(:) = nan;
    
    % Add the confidence interval signs to the source neurons.
    conn_neurons = cellstr(conn_table.Source);
    for i = 1:length(neurons)
        neurons_i = startsWith(conn_neurons, neurons{i});
        conn_table.source_sign(neurons_i) = ci_signs(i);
    end
end

% Group the connectome into a class-only table.
class_table = conn_table(1,[4:5,3,6:end]);
class_table.Properties.VariableNames{1} = 'Source';
class_table.Properties.VariableNames{2} = 'Target';
%class_table = movevars(class_table, 'Weight', 'After', 'Target');
for i = 2:size(conn_table,1)
    
    % Determine the source & target.
    sources = conn_table.source_class(i);
    targets = conn_table.target_class(i);
    
    % Sanity check. 
    row_i = find(class_table.Source == sources & ...
        class_table.Target == targets & ...
        class_table.Type == conn_table.Type(i));
    if length(row_i) > 1
        error('"%s" to "%s" has more than 1 row in the table!', ...
            sources, targets);
    end
    
    % Create the row.
    if isempty(row_i)
        j = size(class_table,1) + 1;
        class_table(j,:) = conn_table(i,[4:5,3,6:end]);
        
    % Add the row info to the existing row.
    else
        % Sum the weights.
        class_table.Weight(row_i) = ...
            class_table.Weight(row_i) + conn_table.Weight(i);
        
        % Sanity check.
        if ~isequaln(class_table(row_i, 4:end), conn_table(i, 6:end))
            error('class_table(%d) ~= conn_table(%d)!', row_i, i);
        end
    end
end

% Sort the tables.
conn_table = sortrows(conn_table, [11, 1:2, 6]);
class_table = sortrows(class_table, [9, 1:2, 4]);
end
