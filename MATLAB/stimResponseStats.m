%% Compute mean sign (stim_sign), t-tests (pT, ci), & wilcoxon ranksum (pW).
% Column = Stimulus
% 1 = butanone
% 2 = pentanedione
% 3 = NaCl
% 4 = light
load('exp_data.mat');

% Which experiments have responses?
is_response = nan(length(exps),3);
for i = 1:length(exps)
    is_response(i,1) = exps(i).response.but;
    is_response(i,2) = exps(i).response.pent;
    is_response(i,3) = exps(i).response.nacl;
    is_response(i,4) = exps(i).response.light;
end
is_response = logical(is_response);

% Compute the neuron stats.
stats.neuron = cell(length(neurons),1);
stats.N = nan(length(neurons),4);
stats.on.mean = nan(length(neurons),4);
stats.off.mean = nan(length(neurons),3);
stats.on.sum_sign = nan(length(neurons),4);
stats.off.sum_sign = nan(length(neurons),3);
stats.on.ci = nan(length(neurons),4,2);
stats.off.ci = nan(length(neurons),3,2);
stats.on.pT = nan(length(neurons),4);
stats.off.pT = nan(length(neurons),3);
stats.on.pW = nan(length(neurons),4);
stats.off.pW = nan(length(neurons),3);
for i = 1:length(neurons)
    
    % Get the neuron name.
    stats.neuron{i} = neurons(i).name;
    
    % Setup the stimulus values.
    if ~neurons(i).isLR
        
        % Signal mean.
        but_pre_mean = neurons(i).but.pre.mean(is_response(:,1));
        but_stim_mean = neurons(i).but.stim.mean(is_response(:,1));
        but_post_mean = neurons(i).but.post.mean(is_response(:,1));
        pent_pre_mean = neurons(i).pent.pre.mean(is_response(:,2));
        pent_stim_mean = neurons(i).pent.stim.mean(is_response(:,2));
        pent_post_mean = neurons(i).pent.post.mean(is_response(:,2));
        nacl_pre_mean = neurons(i).nacl.pre.mean(is_response(:,3));
        nacl_stim_mean = neurons(i).nacl.stim.mean(is_response(:,3));
        nacl_post_mean = neurons(i).nacl.post.mean(is_response(:,3));
        light_stim_mean = neurons(i).light.stim.mean(is_response(:,4));
        light_post_mean = neurons(i).light.post.mean(is_response(:,4));
        
        % Signal start & end means.
        but_stim_end_mean = neurons(i).but.stim.end.mean(is_response(:,1));
        but_post_start_mean = neurons(i).but.post.start.mean(is_response(:,1));
        pent_stim_end_mean = neurons(i).pent.stim.end.mean(is_response(:,2));
        pent_post_start_mean = neurons(i).pent.post.start.mean(is_response(:,2));
        nacl_stim_end_mean = neurons(i).nacl.stim.end.mean(is_response(:,3));
        nacl_post_start_mean = neurons(i).nacl.post.start.mean(is_response(:,3));
        
    % Aggregate the left & right stimulus values.
    else
        
        % Signal mean.
        but_pre_mean = [neurons(i).left.but.pre.mean(is_response(:,1)), ...
            neurons(i).right.but.pre.mean(is_response(:,1))];
        but_stim_mean = [neurons(i).left.but.stim.mean(is_response(:,1)), ...
            neurons(i).right.but.stim.mean(is_response(:,1))];
        but_post_mean = [neurons(i).left.but.post.mean(is_response(:,1)), ...
            neurons(i).right.but.post.mean(is_response(:,1))];
        pent_pre_mean = [neurons(i).left.pent.pre.mean(is_response(:,2)), ...
            neurons(i).right.pent.pre.mean(is_response(:,2))];
        pent_stim_mean = [neurons(i).left.pent.stim.mean(is_response(:,2)), ...
            neurons(i).right.pent.stim.mean(is_response(:,2))];
        pent_post_mean = [neurons(i).left.pent.post.mean(is_response(:,2)), ...
            neurons(i).right.pent.post.mean(is_response(:,2))];
        nacl_pre_mean = [neurons(i).left.nacl.pre.mean(is_response(:,3)), ...
            neurons(i).right.nacl.pre.mean(is_response(:,3))];
        nacl_stim_mean = [neurons(i).left.nacl.stim.mean(is_response(:,3)), ...
            neurons(i).right.nacl.stim.mean(is_response(:,3))];
        nacl_post_mean = [neurons(i).left.nacl.post.mean(is_response(:,3)), ...
            neurons(i).right.nacl.post.mean(is_response(:,3))];
        light_stim_mean = [neurons(i).left.light.stim.mean(is_response(:,4)), ...
            neurons(i).right.light.stim.mean(is_response(:,4))];
        light_post_mean = [neurons(i).left.light.post.mean(is_response(:,4)), ...
            neurons(i).right.light.post.mean(is_response(:,4))];
        
        % Signal start & end means.
        but_stim_end_mean = [neurons(i).left.but.stim.end.mean(is_response(:,1)), ...
            neurons(i).right.but.stim.end.mean(is_response(:,1))];
        but_post_start_mean = [neurons(i).left.but.post.start.mean(is_response(:,1)), ...
            neurons(i).right.but.post.start.mean(is_response(:,1))];
        pent_stim_end_mean = [neurons(i).left.pent.stim.end.mean(is_response(:,2)), ...
            neurons(i).right.pent.stim.end.mean(is_response(:,2))];
        pent_post_start_mean = [neurons(i).left.pent.post.start.mean(is_response(:,2)), ...
            neurons(i).right.pent.post.start.mean(is_response(:,2))];
        nacl_stim_end_mean = [neurons(i).left.nacl.stim.end.mean(is_response(:,3)), ...
            neurons(i).right.nacl.stim.end.mean(is_response(:,3))];
        nacl_post_start_mean = [neurons(i).left.nacl.post.start.mean(is_response(:,3)), ...
            neurons(i).right.nacl.post.start.mean(is_response(:,3))];
    end
    
    % Do we have enough data?
    % Note: testing for the presence of the stim periods is sufficient.
    min_N = 3;
    is_but_mean_data = sum(~isnan(but_pre_mean)) > min_N && ...
        sum(~isnan(but_stim_mean)) > min_N &&  sum(~isnan(but_post_mean)) > min_N;
    is_pent_mean_data = sum(~isnan(pent_pre_mean)) > min_N && ...
        sum(~isnan(pent_stim_mean)) > min_N && sum(~isnan(pent_post_mean)) > min_N;
    is_nacl_mean_data = sum(~isnan(nacl_pre_mean)) > min_N && ...
        sum(~isnan(nacl_stim_mean)) > 2 && sum(~isnan(nacl_post_mean)) > min_N;
    is_light_mean_data = sum(~isnan(light_stim_mean)) > min_N && ...
        sum(~isnan(light_post_mean)) > min_N;
    
    % Compute the mean, sign, & samples.
    if is_but_mean_data
        
        % Butanone on.
        but_on = but_stim_mean - but_pre_mean;
        stats.on.mean(i,1) = nanmean(but_on);
        stats.on.sum_sign(i,1) = nansum(sign(but_on));
        stats.N(i,1) = sum(~isnan(but_stim_mean));

        % Butanone off.
        %but_off = but_post_start_mean - but_stim_end_mean;
        but_off = but_post_mean - but_stim_mean;
        stats.off.mean(i,1) = nanmean(but_off);
        stats.off.sum_sign(i,1) = nansum(sign(but_off));
    end
    if is_pent_mean_data
        
        % Pentanedione on.
        pent_on = pent_stim_mean - pent_pre_mean;
        stats.on.mean(i,2) = nanmean(pent_on);
        stats.on.sum_sign(i,2) = nansum(sign(pent_on));
        stats.N(i,2) = sum(~isnan(pent_stim_mean));
        
        % Pentanedione off.
        %pent_off = pent_post_start_mean - pent_stim_end_mean;
        pent_off = pent_post_mean - pent_stim_mean;
        stats.off.mean(i,2) = nanmean(pent_off);
        stats.off.sum_sign(i,2) = nansum(sign(pent_off));
    end
    if is_nacl_mean_data
        
        % NaCl on.
        nacl_on = nacl_stim_mean - nacl_pre_mean;
        stats.on.mean(i,3) = nanmean(nacl_on);
        stats.on.sum_sign(i,3) = nansum(sign(nacl_on));
        stats.N(i,3) = sum(~isnan(nacl_stim_mean));
        
        % NaCl off.
        %nacl_off = nacl_post_start_mean - nacl_stim_end_mean;
        nacl_off = nacl_post_mean - nacl_stim_mean;
        stats.off.mean(i,3) = nanmean(nacl_off);
        stats.off.sum_sign(i,3) = nansum(sign(nacl_off));
    end
    if is_light_mean_data
        
        % Light on.
        light_on = light_stim_mean - light_post_mean;
        stats.on.mean(i,4) = nanmean(light_on);
        stats.on.sum_sign(i,4) = nansum(sign(light_on));
        stats.N(i,4) = sum(~isnan(light_stim_mean));
    end
    
    % Compute the paired t-tests & confidence intervals.
    if is_but_mean_data
        [~,stats.on.pT(i,1), stats.on.ci(i,1,:)] = ...
            ttest(but_stim_mean, but_pre_mean);
        [~,stats.off.pT(i,1), stats.off.ci(i,1,:)] = ...
            ttest(but_post_mean, but_stim_mean, 'tail', 'right');
        %[~,stats.off.pT(i,1), stats.off.ci(i,1,:)] = ...
            %ttest(but_post_start_mean, but_stim_end_mean);
    end
    if is_pent_mean_data
        [~,stats.on.pT(i,2), stats.on.ci(i,2,:)] = ...
            ttest(pent_stim_mean, pent_pre_mean);
        [~,stats.off.pT(i,2), stats.off.ci(i,2,:)] = ...
            ttest(pent_post_mean, pent_stim_mean, 'tail', 'right');
        %[~,stats.off.pT(i,2), stats.off.ci(i,2,:)] = ...
            %ttest(pent_post_start_mean, pent_stim_end_mean);
    end
    if is_nacl_mean_data
        [~,stats.on.pT(i,3), stats.on.ci(i,3,:)] = ...
            ttest(nacl_stim_mean, nacl_pre_mean);
        [~,stats.off.pT(i,3), stats.off.ci(i,3,:)] = ...
            ttest(nacl_post_mean, nacl_stim_mean, 'tail', 'right');
        %[~,stats.off.pT(i,3), stats.off.ci(i,3,:)] = ...
            %ttest(nacl_post_start_mean, nacl_stim_end_mean);
    end
    if is_light_mean_data
        [~,stats.on.pT(i,4), stats.on.ci(i,4,:)] = ...
            ttest(light_stim_mean, light_post_mean);
    end
    
    % Compute the Wilcoxon rank sums.
    if is_but_mean_data
        stats.on.pW(i,1) = ranksum(but_stim_mean, but_pre_mean);
        stats.off.pW(i,1) = ranksum(but_post_mean, but_stim_mean, ...
            'tail', 'right');
        %stats.off.pW(i,1) = ranksum(but_post_start_mean, but_stim_end_mean);
    end
    if is_pent_mean_data
        stats.on.pW(i,2) = ranksum(pent_stim_mean, pent_pre_mean);
        stats.off.pW(i,2) = ranksum(pent_post_mean, pent_stim_mean, ...
            'tail', 'right');
        %stats.off.pW(i,2) = ranksum(pent_post_start_mean, pent_stim_end_mean);
    end
    if is_nacl_mean_data
        stats.on.pW(i,3) = ranksum(nacl_stim_mean, nacl_pre_mean);
        stats.off.pW(i,3) = ranksum(nacl_post_mean, nacl_stim_mean, ...
            'tail', 'right');
        %stats.off.pW(i,3) = ranksum(nacl_post_start_mean, nacl_stim_end_mean);
    end
    if is_light_mean_data
        stats.on.pW(i,4) = ranksum(light_stim_mean, light_post_mean);
    end
end


%% Compute the neuron signs from the CI (excitation = +1 & inhibition = -1).
stats.on.ci_sign = nan(size(stats.on.ci,1),size(stats.on.ci,2));
stats.off.ci_sign = nan(size(stats.off.ci,1),size(stats.off.ci,2));
for i = 1:size(stats.on.ci,2)
    stats.on.ci_sign((stats.on.ci(:,i,1) >= 0 & stats.on.ci(:,i,2) > 0), i) = 1;
    stats.on.ci_sign((stats.on.ci(:,i,1) < 0 & stats.on.ci(:,i,2) <= 0), i) = -1;
end
for i = 1:size(stats.off.ci,2)
    stats.off.ci_sign((stats.off.ci(:,i,1) >= 0 & stats.off.ci(:,i,2) > 0), i) = 1;
    stats.off.ci_sign((stats.off.ci(:,i,1) < 0 & stats.off.ci(:,i,2) <= 0), i) = -1;
end


%% Compute OFF responses.
% Note: we're comparing 2 transitions: 1) post vs. stim 2) stim vs. pre.
% A. Opposing signs between 1 & 2 is a return to baseline (ignore these).
% B. Inhibition for 1 & 2 means a continuously decreasing signal,
%    could be a long response, overshoot, or bleach trend (ignore these).
% C. Excitation for 1 & 2 means a continuously increasing signal, by eye,
%    these are the only unambiguous looking OFF responses.
keep_i = ((stats.on.mean(:,1:3) > 0) & (stats.off.mean(:,:) > 0));
stats.off.pT(~keep_i) = nan;
stats.off.pW(~keep_i) = nan;


%% Remove the motor neurons/circuit.
load('neuron_props.mat');
table_neurons = cellstr(neuron_table.Neuron);

% Remove 0-padding from the motor neurons.
for i = 1:length(table_neurons)
    name = table_neurons{i};
    if length(name) > 3 && name(3) == '0'
        table_neurons{i} = name([1:2,4:end]);
    end
end

% Find the motor neurons.
is_motor = false(length(stats.neuron),1);
for i = 1:length(stats.neuron)
    
    % Find motor neurons.
    table_i = find(startsWith(table_neurons, stats.neuron{i}),1);
    is_motor(i) = neuron_table(table_i,:).Motor;
    
    % Find command neurons.
    if strcmp(stats.neuron{i}, 'AVA') || ...
            strcmp(stats.neuron{i}, 'AVB') || ...
            strcmp(stats.neuron{i}, 'AVD') || ...
            strcmp(stats.neuron{i}, 'AVE') || ...
            strcmp(stats.neuron{i}, 'PVC')
        is_motor(i) = true;
    end
end

% Remove the motor neurons.
all_stats = stats;
stats.neuron = stats.neuron(~is_motor);
stats.N = stats.N(~is_motor,:);
stats.on.mean = stats.on.mean(~is_motor,:);
stats.off.mean = stats.off.mean(~is_motor,:);
stats.on.sum_sign = stats.on.sum_sign(~is_motor,:);
stats.off.sum_sign = stats.off.sum_sign(~is_motor,:);
stats.on.ci = stats.on.ci(~is_motor,:,:);
stats.off.ci = stats.off.ci(~is_motor,:,:);
stats.on.ci_sign = stats.on.ci_sign(~is_motor,:);
stats.off.ci_sign = stats.off.ci_sign(~is_motor,:);
stats.on.pT = stats.on.pT(~is_motor,:);
stats.off.pT = stats.off.pT(~is_motor,:);
stats.on.pW = stats.on.pW(~is_motor,:);
stats.off.pW = stats.off.pW(~is_motor,:);


%% Correct, within stimulus, for multiple testing.
is_BHFDR = false;

% Stimuli on.
pT_on_data = ~isnan(stats.on.pT);
pW_on_data = ~isnan(stats.on.pW);
stats.on.qT.stim = nan(size(stats.on.pT));
stats.on.qW.stim = nan(size(stats.on.pW));
for i = 1:size(pT_on_data,2)
stats.on.qT.stim(pT_on_data(:,i),i) = ...
    mafdr(stats.on.pT(pT_on_data(:,i),i), 'BHFDR', is_BHFDR);
stats.on.qW.stim(pW_on_data(:,i),i) = ...
    mafdr(stats.on.pW(pW_on_data(:,i),i), 'BHFDR', is_BHFDR);
end

% Stimuli off.
pT_off_data = ~isnan(stats.off.pT);
pW_off_data = ~isnan(stats.off.pW);
stats.off.qT.stim = nan(size(stats.off.pT));
stats.off.qW.stim = nan(size(stats.off.pW));
for i = 1:size(pT_off_data,2)
stats.off.qT.stim(pT_off_data(:,i),i) = ...
    mafdr(stats.off.pT(pT_off_data(:,i),i), 'BHFDR', is_BHFDR);
stats.off.qW.stim(pW_off_data(:,i),i) = ...
    mafdr(stats.off.pW(pW_off_data(:,i),i), 'BHFDR', is_BHFDR);
end


%% Correct, within on/off, for multiple testing.

% Stimuli on.
stats.on.qT.on = nan(size(stats.on.pT));
stats.on.qW.on = nan(size(stats.on.pW));
stats.on.qT.on(pT_on_data) = ...
    mafdr(stats.on.pT(pT_on_data), 'BHFDR', is_BHFDR);
stats.on.qW.on(pW_on_data) = ...
    mafdr(stats.on.pW(pW_on_data), 'BHFDR', is_BHFDR);

% Stimuli off.
stats.off.qT.off = nan(size(stats.off.pT));
stats.off.qW.off = nan(size(stats.off.pW));
stats.off.qT.off(pT_off_data) = ...
    mafdr(stats.off.pT(pT_off_data), 'BHFDR', is_BHFDR);
stats.off.qW.off(pW_off_data) = ...
    mafdr(stats.off.pW(pW_off_data), 'BHFDR', is_BHFDR);


%% Correct across all multiple testing.
stats.on.qT.all = nan(size(stats.on.pT));
stats.off.qT.all = nan(size(stats.off.pT));
stats.on.qW.all = nan(size(stats.on.pW));
stats.off.qW.all = nan(size(stats.off.pW));
all_pT(:,1:4) = stats.on.pT;
all_pT(:,5:7) = stats.off.pT;
all_pW(:,1:4) = stats.on.pW;
all_pW(:,5:7) = stats.off.pW;
pT_data = ~isnan(all_pT);
pW_data = ~isnan(all_pW);
all_pT(pT_data) = mafdr(all_pT(pT_data), 'BHFDR', is_BHFDR);
all_pW(pW_data) = mafdr(all_pW(pW_data), 'BHFDR', is_BHFDR);
stats.on.qT.all(:,1:4) = all_pT(:,1:4);
stats.off.qT.all(:,1:3) = all_pT(:,5:7);
stats.on.qW.all(:,1:4) = all_pW(:,1:4);
stats.off.qW.all(:,1:3) = all_pW(:,5:7);


%% Save the data.
save('exp_stats.mat', 'stats', 'all_stats');


%% Extract the connectome for each stimulus.
load('connectome_syn_gj.mat', 'conn_edges');
sources = cellstr(conn_edges.Source);
targets = cellstr(conn_edges.Target);
premotors = categorical({'AVA', 'AVB', 'AVD', 'AVE', 'PVC'})';

% q-value significance threshold.
q_thresh = 0.1;

% Extract the premotor targets & their motor targets.
premotor_source_i = [];
premotor_target_i = [];
premotor2motor_i = [];
motor2premotor_i = [];
for i = 1:length(premotors)
    premotor_source_i = [premotor_source_i; ...
        find(conn_edges.source_class == premotors(i))];
    premotor_target_i = [premotor_target_i; ...
        find(conn_edges.target_class == premotors(i))];
    premotor2motor_i = [premotor2motor_i; ...
        find(conn_edges.source_class == premotors(i) & conn_edges.target_motor)];
    motor2premotor_i = [motor2premotor_i; ...
        find(conn_edges.target_class == premotors(i) & conn_edges.source_motor)];
end
premotor_source_i = unique(premotor_source_i);
premotor_target_i = unique(premotor_target_i);
connectome.motor.premotor.source.i = premotor_source_i;
connectome.motor.premotor.target.i = premotor_source_i;
connectome.motor.premotor.premotor2motor.i = premotor2motor_i;
connectome.motor.premotor.motor2premotor.i = premotor2motor_i;


%% Extract the butanone response connectome.

% Butanone ON connectome.
but_on_q_i = stats.on.qT.on(:,1) < q_thresh;
but_on_neurons = stats.neuron(but_on_q_i);
connectome.but.on.neurons = but_on_neurons;
but_on_neurons = [but_on_neurons; cellstr(premotors)];
but_on_source = [];
but_on_target = [];
for i = 1:length(but_on_neurons)
    but_on_source = ...
        [but_on_source; find(startsWith(sources, but_on_neurons{i}))];
    but_on_target = ...
        [but_on_target; find(startsWith(targets, but_on_neurons{i}))];
end
but_on_conn_i = intersect(but_on_source, but_on_target);
connectome.but.on.i = but_on_conn_i;

% Sanity check the butanone ON response connectome.
but_on_conn_neurons = cellstr(union(conn_edges.Source(but_on_conn_i), ...
    conn_edges.Target(but_on_conn_i)));
for i = 1:length(but_on_neurons)
    if ~any(startsWith(but_on_conn_neurons, but_on_neurons{i}))
        warning('"%s" is not represented in the butanone ON connectome!', ...
            but_on_neurons{i});
    end
end

% Butanone OFF connectome.
but_off_q_i = stats.off.qT.off(:,1) < q_thresh;
but_off_neurons = stats.neuron(but_off_q_i);
connectome.but.off.neurons = but_off_neurons;
but_off_neurons = [but_off_neurons; cellstr(premotors)];
but_off_source = [];
but_off_target = [];
for i = 1:length(but_off_neurons)
    but_off_source = ...
        [but_off_source; find(startsWith(sources, but_off_neurons{i}))];
    but_off_target = ...
        [but_off_target; find(startsWith(targets, but_off_neurons{i}))];
end
but_off_conn_i = intersect(but_off_source, but_off_target);
connectome.but.off.i = but_off_conn_i;

% Sanity check the butanone OFF response connectome.
but_off_conn_neurons = cellstr(union(conn_edges.Source(but_off_conn_i), ...
    conn_edges.Target(but_off_conn_i)));
for i = 1:length(but_off_neurons)
    if ~any(startsWith(but_off_conn_neurons, but_off_neurons{i}))
        warning('"%s" is not represented in the butanone OFF connectome!', ...
            but_off_neurons{i});
    end
end


%% Extract the pentanedione response connectome.

% Pentanedione ON connectome.
pent_on_q_i = stats.on.qT.on(:,2) < q_thresh;
pent_on_neurons = stats.neuron(pent_on_q_i);
connectome.pent.on.neurons = pent_on_neurons;
pent_on_neurons = [pent_on_neurons; cellstr(premotors)];
pent_on_source = [];
pent_on_target = [];
for i = 1:length(pent_on_neurons)
    pent_on_source = ...
        [pent_on_source; find(startsWith(sources, pent_on_neurons{i}))];
    pent_on_target = ...
        [pent_on_target; find(startsWith(targets, pent_on_neurons{i}))];
end
pent_on_conn_i = intersect(pent_on_source, pent_on_target);
connectome.pent.on.i = pent_on_conn_i;

% Sanity check the pentanedione ON response connectome.
pent_on_conn_neurons = cellstr(union(conn_edges.Source(pent_on_conn_i), ...
    conn_edges.Target(pent_on_conn_i)));
for i = 1:length(pent_on_neurons)
    if ~any(startsWith(pent_on_conn_neurons, pent_on_neurons{i}))
        warning('"%s" is not represented in the pentanedione ON connectome!', ...
            pent_on_neurons{i});
    end
end

% Pentanedione OFF connectome.
pent_off_q_i = stats.off.qT.off(:,2) < q_thresh;
pent_off_neurons = stats.neuron(pent_off_q_i);
connectome.pent.off.neurons = pent_off_neurons;
pent_off_neurons = [pent_off_neurons; cellstr(premotors)];
pent_off_source = [];
pent_off_target = [];
for i = 1:length(pent_off_neurons)
    pent_off_source = ...
        [pent_off_source; find(startsWith(sources, pent_off_neurons{i}))];
    pent_off_target = ...
        [pent_off_target; find(startsWith(targets, pent_off_neurons{i}))];
end
pent_off_conn_i = intersect(pent_off_source, pent_off_target);
connectome.pent.off.i = pent_off_conn_i;

% Sanity check the pentanedione OFF response connectome.
pent_off_conn_neurons = cellstr(union(conn_edges.Source(pent_off_conn_i), ...
    conn_edges.Target(pent_off_conn_i)));
for i = 1:length(pent_off_neurons)
    if ~any(startsWith(pent_off_conn_neurons, pent_off_neurons{i}))
        warning('"%s" is not represented in the pentanedione OFF connectome!', ...
            pent_off_neurons{i});
    end
end


%% Extract the NaCl response connectome.

% NaCl ON connectome.
nacl_on_q_i = stats.on.qT.on(:,3) < q_thresh;
nacl_on_neurons = stats.neuron(nacl_on_q_i);
connectome.nacl.on.neurons = nacl_on_neurons;
nacl_on_neurons = [nacl_on_neurons; cellstr(premotors)];
nacl_on_source = [];
nacl_on_target = [];
for i = 1:length(nacl_on_neurons)
    nacl_on_source = ...
        [nacl_on_source; find(startsWith(sources, nacl_on_neurons{i}))];
    nacl_on_target = ...
        [nacl_on_target; find(startsWith(targets, nacl_on_neurons{i}))];
end
nacl_on_conn_i = intersect(nacl_on_source, nacl_on_target);
connectome.nacl.on.i = nacl_on_conn_i;

% Sanity check the NaCl ON response connectome.
nacl_on_conn_neurons = cellstr(union(conn_edges.Source(nacl_on_conn_i), ...
    conn_edges.Target(nacl_on_conn_i)));
for i = 1:length(nacl_on_neurons)
    if ~any(startsWith(nacl_on_conn_neurons, nacl_on_neurons{i}))
        warning('"%s" is not represented in the NaCl ON connectome!', ...
            nacl_on_neurons{i});
    end
end

% NaCl OFF connectome.
nacl_off_q_i = stats.off.qT.off(:,3) < q_thresh;
nacl_off_neurons = stats.neuron(nacl_off_q_i);
connectome.nacl.off.neurons = nacl_off_neurons;
nacl_off_neurons = [nacl_off_neurons; cellstr(premotors)];
nacl_off_source = [];
nacl_off_target = [];
for i = 1:length(nacl_off_neurons)
    nacl_off_source = ...
        [nacl_off_source; find(startsWith(sources, nacl_off_neurons{i}))];
    nacl_off_target = ...
        [nacl_off_target; find(startsWith(targets, nacl_off_neurons{i}))];
end
nacl_off_conn_i = intersect(nacl_off_source, nacl_off_target);
connectome.nacl.off.i = nacl_off_conn_i;

% Sanity check the NaCl OFF response connectome.
nacl_off_conn_neurons = cellstr(union(conn_edges.Source(nacl_off_conn_i), ...
    conn_edges.Target(nacl_off_conn_i)));
for i = 1:length(nacl_off_neurons)
    if ~any(startsWith(nacl_off_conn_neurons, nacl_off_neurons{i}))
        warning('"%s" is not represented in the NaCl OFF connectome!', ...
            nacl_off_neurons{i});
    end
end


%% Extract the light response connectome.
light_on_q_i = stats.on.qT.on(:,4) < q_thresh;
light_on_neurons = stats.neuron(light_on_q_i);
connectome.light.on.neurons = light_on_neurons;
light_on_neurons = [light_on_neurons; cellstr(premotors)];
light_on_source = [];
light_on_target = [];
for i = 1:length(light_on_neurons)
    light_on_source = ...
        [light_on_source; find(startsWith(sources, light_on_neurons{i}))];
    light_on_target = ...
        [light_on_target; find(startsWith(targets, light_on_neurons{i}))];
end
light_on_conn_i = intersect(light_on_source, light_on_target);
connectome.light.on.i = light_on_conn_i;

% Sanity check the light response connectome.
light_on_conn_neurons = cellstr(union(conn_edges.Source(light_on_conn_i), ...
    conn_edges.Target(light_on_conn_i)));
for i = 1:length(light_on_neurons)
    if ~any(startsWith(light_on_conn_neurons, light_on_neurons{i}))
        warning('"%s" is not represented in the light ON connectome!', ...
            light_on_neurons{i});
    end
end


%% Construct the connectomes.

% Construct the ON connectomes.
[connectome.but.on.edges, connectome.but.on.class_edges] = ...
    buildConnTable(conn_edges, connectome.but.on.i, ...
    stats.neuron(but_on_q_i), stats.on.ci_sign(but_on_q_i,1), true);
[connectome.pent.on.edges, connectome.pent.on.class_edges] = ...
    buildConnTable(conn_edges, connectome.pent.on.i, ...
    stats.neuron(pent_on_q_i), stats.on.ci_sign(pent_on_q_i,2), true);
[connectome.nacl.on.edges, connectome.nacl.on.class_edges] = ...
    buildConnTable(conn_edges, connectome.nacl.on.i, ...
    stats.neuron(nacl_on_q_i), stats.on.ci_sign(nacl_on_q_i,3), true);
[connectome.light.on.edges, connectome.light.on.class_edges] = ...
    buildConnTable(conn_edges, connectome.light.on.i, ...
    stats.neuron(light_on_q_i), stats.on.ci_sign(light_on_q_i,4), true);

% Construct the OFF connectomes.
[connectome.but.off.edges, connectome.but.off.class_edges] = ...
    buildConnTable(conn_edges, connectome.but.off.i, ...
    stats.neuron(but_off_q_i), stats.off.ci_sign(but_off_q_i,1), true);
[connectome.pent.off.edges, connectome.pent.off.class_edges] = ...
    buildConnTable(conn_edges, connectome.pent.off.i, ...
    stats.neuron(pent_off_q_i), stats.off.ci_sign(pent_off_q_i,2), true);
[connectome.nacl.off.edges, connectome.nacl.off.class_edges] = ...
    buildConnTable(conn_edges, connectome.nacl.off.i, ...
    stats.neuron(nacl_off_q_i), stats.off.ci_sign(nacl_off_q_i,3), true);

% Save the ON connectomes for Cytoscape.
writetable(connectome.but.on.edges, 'butanone_ON.csv');
writetable(connectome.but.on.class_edges, 'butanone_class_ON.csv');
writetable(connectome.pent.on.edges, 'pentanedione_ON.csv');
writetable(connectome.pent.on.class_edges, 'pentanedione_class_ON.csv');
writetable(connectome.nacl.on.edges, 'nacl_ON.csv');
writetable(connectome.nacl.on.class_edges, 'nacl_class_ON.csv');
writetable(connectome.light.on.edges, 'light_ON.csv');
writetable(connectome.light.on.class_edges, 'light_class_ON.csv');

% Save the OFF connectomes for Cytoscape.
writetable(connectome.but.off.edges, 'butanone_OFF.csv');
writetable(connectome.but.off.class_edges, 'butanone_class_OFF.csv');
writetable(connectome.pent.off.edges, 'pentanedione_OFF.csv');
writetable(connectome.pent.off.class_edges, 'pentanedione_class_OFF.csv');
writetable(connectome.nacl.off.edges, 'nacl_OFF.csv');
writetable(connectome.nacl.off.class_edges, 'nacl_class_OFF.csv');

% Save the connectomes.
save('exp_conns.mat', 'connectome', 'conn_edges');
