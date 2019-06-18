function [conn_dist, syn_dist, gj_dist] = connDist(conn, neurons)
%CONNDIST The connectome distance (hops) from a set of start neurons.
%
%   Input:
%       conn    - the connectome
%       neurons - the start neurons
%
%   Output:
%       conn_dist - the cumulative neurons at each hop from the start
%       syn_dist - the cumulative neurons at each chemical hop from the start
%       gj_dist - the cumulative neurons at each electrical hop from the start

% Convert the neurons to categoricals.
if ~iscategorical(neurons)
    neurons = categorical(cellstr(neurons));
end

% Traverse the connectome.
is_more = true;
conn_dist.neuron{1} = unique(neurons);
conn_dist.class = [];
while is_more
    
    % Compute the start.
    start_neurons = conn_dist.neuron{end};
    start_i = conn.Source == start_neurons(1);
    for i = 2:length(start_neurons)
        start_i = start_i | conn.Source == start_neurons(i);
    end

    % Add the source classes.
    % Note: categoricals aren't sorted alphabetically.
    conn_classes = unique(conn.source_class(start_i));
    conn_classes(isundefined(conn_classes)) = [];
    conn_classes = categorical(sort(cellstr(conn_classes)));
    conn_dist.class{end + 1} = conn_classes;
    
    % Add the connected neurons.
    % Note: categoricals aren't sorted alphabetically.
    conn_neurons = unique(conn.Target(start_i & conn.target_neuron));
    conn_neurons = union(start_neurons, conn_neurons);
    conn_neurons = categorical(sort(cellstr(conn_neurons)));
    
    % Did we add neurons?
    is_more = length(start_neurons) < length(conn_neurons);
    if is_more
        conn_dist.neuron{end + 1} = conn_neurons;
    end
end

% Traverse the chemical connectome.
is_more = true;
syn_dist.neuron{1} = unique(neurons);
syn_dist.class = [];
while is_more
    
    % Compute the start.
    start_neurons = syn_dist.neuron{end};
    start_i = conn.Source == start_neurons(1);
    for i = 2:length(start_neurons)
        start_i = start_i | conn.Source == start_neurons(i);
    end

    % Add the source classes.
    % Note: categoricals aren't sorted alphabetically.
    conn_classes = unique(conn.source_class(start_i));
    conn_classes(isundefined(conn_classes)) = [];
    conn_classes = categorical(sort(cellstr(conn_classes)));
    syn_dist.class{end + 1} = conn_classes;
    
    % Add the connected neurons.
    % Note: categoricals aren't sorted alphabetically.
    conn_neurons = unique(conn.Target(start_i & conn.target_neuron & ...
        conn.Type == 'chemical'));
    conn_neurons = union(start_neurons, conn_neurons);
    conn_neurons = categorical(sort(cellstr(conn_neurons)));
    
    % Did we add neurons?
    is_more = length(start_neurons) < length(conn_neurons);
    if is_more
        syn_dist.neuron{end + 1} = conn_neurons;
    end
end

% Traverse the chemical connectome.
is_more = true;
gj_dist.neuron{1} = unique(neurons);
gj_dist.class = [];
while is_more
    
    % Compute the start.
    start_neurons = gj_dist.neuron{end};
    start_i = conn.Source == start_neurons(1);
    for i = 2:length(start_neurons)
        start_i = start_i | conn.Source == start_neurons(i);
    end

    % Add the source classes.
    % Note: categoricals aren't sorted alphabetically.
    conn_classes = unique(conn.source_class(start_i));
    conn_classes(isundefined(conn_classes)) = [];
    conn_classes = categorical(sort(cellstr(conn_classes)));
    gj_dist.class{end + 1} = conn_classes;
    
    % Add the connected neurons.
    % Note: categoricals aren't sorted alphabetically.
    conn_neurons = unique(conn.Target(start_i & conn.target_neuron & ...
        conn.Type == 'electrical'));
    conn_neurons = union(start_neurons, conn_neurons);
    conn_neurons = categorical(sort(cellstr(conn_neurons)));
    
    % Did we add neurons?
    is_more = length(start_neurons) < length(conn_neurons);
    if is_more
        gj_dist.neuron{end + 1} = conn_neurons;
    end
end
end
