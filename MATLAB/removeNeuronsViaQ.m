% Setup the trace file.
file = '66_head_run101_traces.mat';

% List of neurons to remove via Q-ing.
q = {
'CEPDL'
'CEPDR'
'IL2R'
'IL1VR'
'SMBDR'
'RIR'
'ADFL'
'AINR'
'ASGL'
'ASJL'
'ASKR'
'AVBL'
'AVDL'
'AVEL'
'RIAR'
'RIBR'
'RIMR'
    };

% Load the neuron names.
load(file, 'neuron_names');
neuron_names = cellstr(neuron_names);

% Q-out the bad neurons.
disp(['Removing ' num2str(length(q)) ' neurons.']);
for i = 1:length(q)
    
    % Find the neuron.
   neuron_i = find(startsWith(cellstr(neuron_names), q{i}));
   if length(neuron_i) ~= 1
       error('Cannot find "%s"!', q{i});
   end
   
   % Q-out the neuron.
   name = ['Q_' neuron_names{neuron_i}];
   disp(['Renaming "' q{i} '" to "' name '".']);
    neuron_names{neuron_i} = name;
end

% Save the neurons to the file.
save(file, 'neuron_names', '-append');
