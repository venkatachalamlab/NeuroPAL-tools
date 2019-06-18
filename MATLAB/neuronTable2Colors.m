function neurons = neuronTable2Colors(neuron_table)
%NEURONTABLE2COLORS Convert a neuron table to a list of neurons & their colors.
%   [NEURONS COLORS] = NEURONTABLE2COLORS(NEURON_TABLE)
%
%   Input:
%   neuron_table - a table of neurons & their colors
%
%   Outputs:
%   neurons - a structure of neurons & their colors


% Initilaize the table offsets.
start_off = 2;
end_off = 1;

% Convert the table to a struct list.
num_neurons = size(neuron_table,1) - (start_off + end_off);
neurons(1:num_neurons) = struct('name', string([]), 'RGB', zeros(1,3));
for i = 1:num_neurons
    iTable = i + start_off;
    neurons(i).name = neuron_table.Neurons{iTable};
    neurons(i).RGB(1) = neuron_table.RedValue(iTable)/255;
    neurons(i).RGB(2) = neuron_table.GreenValue(iTable)/255;
    neurons(i).RGB(3) = neuron_table.BlueValue(iTable)/255;
end

% Get rid of duplicates.
[~,i,~] = unique({neurons.name});
neurons = neurons(i);

% Expand neuron lists.
old_neurons = neurons;
iNeurons = 1;
for i = 1:length(old_neurons)
    
    % Do we have  "-" or "/"?
    iDash = strfind(old_neurons(i).name, '-');
    iSlash = strfind(old_neurons(i).name, '/');
    if isempty(iDash) && isempty(iSlash)
        neurons(iNeurons) = old_neurons(i);
        iNeurons = iNeurons + 1;
    
    % Separate the neurons.
    elseif isempty(iDash)
        
        % Neuron before slash.
        old_name = old_neurons(i).name;
        neurons(iNeurons).name = old_name(1:(iSlash - 1));
        neurons(iNeurons).RGB = old_neurons(i).RGB;
        iNeurons = iNeurons + 1;
        
        % Neuron after slash.
        neurons(iNeurons).name = [old_name(1:(iSlash - 2))...
            old_name((iSlash + 1):end)];
        neurons(iNeurons).RGB = old_neurons(i).RGB;
        iNeurons = iNeurons + 1;
        
    % Expand the neurons list.
    else
        
        % Find the numbers in the neuron name.
        old_name = old_neurons(i).name;
        rangeStart = iDash - 1;
        while rangeStart > 1 && isstrprop(old_name(rangeStart), 'digit')
            rangeStart = rangeStart - 1;
        end
        range = sscanf(old_name((rangeStart+1):end), '%d-%d');
        
        % Expand the neuron list.
        name = old_name(1:(rangeStart));
        for j = range(1):range(2)
            neurons(iNeurons).name = [name num2str(j)];
            neurons(iNeurons).RGB = old_neurons(i).RGB;
            iNeurons = iNeurons + 1;
        end
    end
end

% Sort the neurons.
[~, iNeurons] = natsort({neurons.name});
neurons = neurons(iNeurons);
end
