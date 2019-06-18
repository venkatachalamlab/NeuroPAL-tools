%% Build the connectome from the chemical synapses & gap juctions.

% Load the the connectomes.
%chem = readtable('herm_chem.csv');
load('Metabotropic_NTRs.mat', 'mNTR_edges');
gj = readtable('herm_gj.csv');
gj.Properties.VariableNames{1} = 'Source';

% Clean up the sources.
gj.Source = categorical(gj.Source);

% The existing connections are chemical (we need to add the gap junctions).
mNTR_edges.Type(:) = categorical(cellstr('chemical'));
conn_edges = mNTR_edges;

% Convert the gap juction matrix to a table. 
electrical = categorical(cellstr('electrical'));
for i=1:size(gj,1)
    for j = 2:size(gj,2)
        
        % Convert the weight to a number.
        weight = gj{i,j};
        if isstring(weight)
            weight = str2double(weight);
        elseif iscategorical(weight)
            weight = str2double(cellstr(weight));
        end
        
        % Is the weight positive?
        if isnan(weight) || ~isnumeric(weight) || weight < 1
            weight = [];
        end
        
        % Does the edge have a weight?
        if ~isempty(weight)
            
            % Did we already add the GJ edge?
            % Note: I want to keep the table compact, one entry per edge.
            % GJ edges are non-directional with single row representation.
            source = gj.Source(i);
            target = categorical(cellstr(gj.Properties.VariableNames{j}));
            if any((conn_edges.Source == target & ...
                    conn_edges.Target == source) & ...
                    conn_edges.Type == electrical)
                continue;
            end
            
            % Find the source.
            source_i = find(source == mNTR_edges.Source, 1);
            source_row = mNTR_edges(source_i,:);
            
            % Is there any info on the source?
            k = size(conn_edges,1) + 1;
            source_neuron = false;
            if isempty(source_i)
                conn_edges.Source(k) = source;
                conn_edges.source_neuron(k) = false;
                conn_edges.neuron_edge(k) = false;
                conn_edges.source_class(k) = categorical(nan);
                conn_edges.source_sensory(k) = false;
                conn_edges.source_inter(k) = false;
                conn_edges.source_motor(k) = false;
                conn_edges.source_poly(k) = false;
                conn_edges.source_pharynx(k) = false;
                conn_edges.classical_NT(k) = false;
                conn_edges.monoamine_NT(k) = false;
                conn_edges.multi_NT(k) = false;
                conn_edges.acetylcholine(k) = false;
                conn_edges.GABA(k) = false;
                conn_edges.glutamate(k) = false;
                conn_edges.dopamine(k) = false;
                conn_edges.octopamine(k) = false;
                conn_edges.serotonin(k) = false;
                conn_edges.tyramine(k) = false;
                conn_edges.orphan(k) = false;
                conn_edges.source_mNTR(k) = false;
                conn_edges.source_GAR(k) = false;
                conn_edges.source_GBB(k) = false;
                conn_edges.source_MGL(k) = false;
                conn_edges.source_GABAa(k) = false;
                
            % Add the source info.
            else
                conn_edges(k,:) = source_row;
                source_neuron = source_row.source_neuron;
            end
            
            % Find the target.
            % Note: all neurons are sources but not all are targets.
            % *** Except CAN.
            target_i = find(target == mNTR_edges.Target, 1);
            target_row = mNTR_edges(target_i,:);
            
            % Add the row.
            conn_edges.Target(k) = target;
            conn_edges.Weight(k) = weight;
            conn_edges.Type(k) = electrical;
            conn_edges.mNTR_edge(k) = false;
            conn_edges.auto_edge(k) = (source == target);
            
            % We have no target info.
            % Note: all neurons are sources but not all are targets.
            % *** Except CAN.
            if isempty(target_i)
                
                % Use the source info.
                target_i = find(target == mNTR_edges.Source, 1);
                
                % We have no source info for the target.
                if isempty(target_i)
                    conn_edges.target_neuron(k) = false;
                    conn_edges.neuron_edge(k) = false;
                    conn_edges.target_class(k) = categorical(nan);
                    conn_edges.target_sensory(k) = false;
                    conn_edges.target_inter(k) = false;
                    conn_edges.target_motor(k) = false;
                    conn_edges.target_poly(k) = false;
                    conn_edges.target_pharynx(k) = false;
                    conn_edges.target_mNTR(k) = false;
                    conn_edges.target_GAR(k) = false;
                    conn_edges.target_GBB(k) = false;
                    conn_edges.target_MGL(k) = false;
                    conn_edges.target_GABAa(k) = false;
                    
                % We can use the source info for the target.
                else
                    target_row = mNTR_edges(target_i,:);
                    target_neuron = target_row.source_neuron;
                    conn_edges.target_neuron(k) = target_neuron;
                    conn_edges.neuron_edge(k) = source_neuron && target_neuron;
                    conn_edges.target_class(k) = target_row.source_class;
                    conn_edges.target_sensory(k) = target_row.source_sensory;
                    conn_edges.target_inter(k) = target_row.source_inter;
                    conn_edges.target_motor(k) = target_row.source_motor;
                    conn_edges.target_poly(k) = target_row.source_poly;
                    conn_edges.target_pharynx(k) = target_row.source_pharynx;
                    conn_edges.target_mNTR(k) = target_row.source_mNTR;
                    conn_edges.target_GAR(k) = target_row.source_GAR;
                    conn_edges.target_GBB(k) = target_row.source_GBB;
                    conn_edges.target_MGL(k) = target_row.source_MGL;
                    conn_edges.target_GABAa(k) = target_row.source_GABAa;
                end
                
                % Initialize the remaining target info.
                conn_edges.target_GABAa_exp1(k) = false;
                conn_edges.target_GABAa_gab1(k) = false;
                conn_edges.target_GABAa_lgc35(k) = false;
                conn_edges.target_GABAa_lgc36(k) = false;
                conn_edges.target_GABAa_lgc37(k) = false;
                conn_edges.target_GABAa_lgc38(k) = false;
                conn_edges.target_GABAa_unc49(k) = false;
            
            % We have the target info.
            else
                target_neuron = target_row.target_neuron;
                conn_edges.target_neuron(k) = target_neuron;
                conn_edges.neuron_edge(k) = source_neuron && target_neuron;
                conn_edges.target_class(k) = target_row.target_class;
                conn_edges.target_sensory(k) = target_row.target_sensory;
                conn_edges.target_inter(k) = target_row.target_inter;
                conn_edges.target_motor(k) = target_row.target_motor;
                conn_edges.target_poly(k) = target_row.target_poly;
                conn_edges.target_pharynx(k) = target_row.target_pharynx;
                conn_edges.target_mNTR(k) = target_row.target_mNTR;
                conn_edges.target_GAR(k) = target_row.target_GAR;
                conn_edges.target_GBB(k) = target_row.target_GBB;
                conn_edges.target_MGL(k) = target_row.target_MGL;
                conn_edges.target_GABAa(k) = target_row.target_GABAa;
                conn_edges.target_GABAa_exp1(k) = target_row.target_GABAa_exp1;
                conn_edges.target_GABAa_gab1(k) = target_row.target_GABAa_gab1;
                conn_edges.target_GABAa_lgc35(k) = target_row.target_GABAa_lgc35;
                conn_edges.target_GABAa_lgc36(k) = target_row.target_GABAa_lgc36;
                conn_edges.target_GABAa_lgc37(k) = target_row.target_GABAa_lgc37;
                conn_edges.target_GABAa_lgc38(k) = target_row.target_GABAa_lgc38;
                conn_edges.target_GABAa_unc49(k) = target_row.target_GABAa_unc49;
            end
        end
    end
end

% Save the connectome.
save('connectome_syn_gj.mat', 'conn_edges');
